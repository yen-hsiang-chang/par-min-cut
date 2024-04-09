#include <iostream>
#include <random>
#include <omp.h>

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/random.h"
#include "gbbs/gbbs/gbbs.h"
#include "gbbs/gbbs/io.h"
#include "gbbs/macros.h"
#include "gbbs/benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "rctree/contraction/contraction.h"
#include "utils/utils.h"

int main(int argc, char* argv[]) {
  // macros
  using W = unsigned int;

  // random number generator
  const size_t seed = 42;
  parlay::random_generator gen(seed);

  // read the graph, assume {u, v} only has (u, v) or (v, u)
  auto G_dir = gbbs::gbbs_io::read_weighted_symmetric_graph<W>(argv[1], false, false);

  // Construct the undirected version
  auto G_edge_list = gbbs::to_edge_array<W>(G_dir);
  auto G_undir = gbbs::gbbs_io::edge_list_to_symmetric_graph<W>(G_edge_list);

  // Feed exponential distribution to MST solver with (v, Exp(w))
  auto E_MST_input = parlay::sequence<std::tuple<uintV, double>>(G_dir.m);
  auto V_MST_input = parlay::sequence<gbbs::vertex_data>(G_dir.n);
  auto G_MST_input = gbbs::symmetric_graph<gbbs::symmetric_vertex, double>(
    V_MST_input.data(), G_dir.n, G_dir.m, [](){}, E_MST_input.data());
  
  // Construct vertex weight for Contraction RCTree
  auto contraction_vertex_weight = parlay::sequence<Contraction_Type<W>>(G_undir.n);
  auto G_weighted_degree = parlay::sequence<W>(G_undir.n);
  parlay::parallel_for(0, G_undir.n, [&](uintV i) {
    auto get_weight = [&](uintV j) {
      return std::get<1>(G_undir.e0[G_undir.v_data[i].offset + j]); 
    };
    auto weights = parlay::delayed_seq<W>(G_undir.v_data[i].degree, get_weight);
    G_weighted_degree[i] = parlay::reduce(weights);
    contraction_vertex_weight[i] = Contraction_Type<W>(G_weighted_degree[i], 0, 0, 0, true);
  });

  // Construct edge weight for Contraction RCTree
  auto contraction_edge_weight = parlay::sequence<Contraction_Type<W>>::from_function(G_undir.n - 1, [&](uintV i) {
    return Contraction_Type<W>(0, 0, i, i, false);
  });

  // Base case where each vertex is a partition
  W G_min_weighted_degree = parlay::reduce(G_weighted_degree, parlay::minimum<W>());

  std::cout << G_min_weighted_degree << std::endl;

  // TODO: determine the number of iterations
  for(int iter = 0; iter < 1; ++iter) {
    // Step 1: Generate weighted random edge ordering via exponential distribution
    auto weighted_ordering_time = -omp_get_wtime();
    gen(); // Advance one step for the rng
    parlay::parallel_for(0, G_dir.m, [&](uintE i) {
      auto [v, w] = G_dir.e0[i];
      std::exponential_distribution<double> exponential(w);
      auto local_rand = gen[i];
      E_MST_input[i] = std::make_tuple(v, exponential(local_rand));
    });
    weighted_ordering_time += omp_get_wtime();

    std::cout << "Weighted Random Edge Ordering Time: "
              << weighted_ordering_time << " seconds" << std::endl;
  
    // Step 2: Compute MST based on weighted random edge ordering
    auto mst_time = -omp_get_wtime();
    // These will be overwritten in MST algo, reset every iteration
    G_MST_input.n = G_dir.n;
    G_MST_input.m = G_dir.m;
    parlay::parallel_for(0, G_dir.n, [&](uintV i) {
      V_MST_input[i] =  G_dir.v_data[i];
    });
    auto E_MST = gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G_MST_input);
    mst_time += omp_get_wtime();

    std::cout << "MST Time: "
              << mst_time << " seconds" << std::endl;
    
    // If it is a forest, the minimum cut is just 0.
    if (E_MST.size() + 1 != G_dir.n) {
      std::cout << "min cut = 0." << std::endl; 
      break;
    }

    // Step 3: Build contraction RCTree
    auto contraction_rctree_time = -omp_get_wtime();
    utils::general_sort(E_MST, [&](const std::tuple<uintV, uintV, double>& lhs, const std::tuple<uintV, uintV, double>& rhs) {
      return std::get<2>(lhs) < std::get<2>(rhs);
    });

    auto MST_edge_list = parlay::sequence<std::pair<uintV, uintV>>::from_function(E_MST.size(), [&](uintV i) {
      auto [u, v, _] = E_MST[i];
      return std::make_pair(u, v);
    });

    auto contraction_rctree = Contraction_RCTree<W>(G_undir.n, MST_edge_list, gen, 
                                                    contraction_vertex_weight, contraction_edge_weight);
    contraction_rctree_time += omp_get_wtime();

    std::cout << "Contraction RCTree Time: "
              << contraction_rctree_time << " seconds" << std::endl;

    // Step 4: Batch operations
    auto contraction_batched_ops_time = -omp_get_wtime();
    auto contraction_mixop = parlay::sequence<Contraction_MixOp<W>>((G_edge_list.size() + MST_edge_list.size() - 1) * 2);
    // Subtract vertex weight
    parlay::parallel_for(0, G_edge_list.E.size(), [&](uintE i) {
      auto [u, v, w] = G_edge_list.E[i];
      contraction_mixop[i * 2] = Contraction_MixOp<W>(contraction_rctree.query_max_path(u, v) * 4 + 1, 0, u, w);
      contraction_mixop[i * 2 + 1] = Contraction_MixOp<W>(contraction_rctree.query_max_path(u, v) * 4 + 2, 0, v, w);
    });
    // Join edges and Query component induced by joined edges
    parlay::parallel_for(0, MST_edge_list.size() - 1, [&](uintE i) {
      auto [u, v] = MST_edge_list[i];
      contraction_mixop[G_edge_list.size() * 2 + i * 2] = Contraction_MixOp<W>(i * 4 + 3, 1, i, 0);
      contraction_mixop[G_edge_list.size() * 2 + i * 2 + 1] = Contraction_MixOp<W>(i * 4 + 4, 2, u, 0);
    });
    // W answer = std::min(G_min_weighted_degree, contraction_rctree.batch_operations_sequential(contraction_mixop));
    // W answer = std::min(G_min_weighted_degree, contraction_rctree.batch_operations(contraction_mixop));
    W answer = contraction_rctree.batch_operations(contraction_mixop);
    contraction_batched_ops_time += omp_get_wtime();

    std::cout << "Contraction Batched Operations Time: "
              << contraction_batched_ops_time << " seconds" << std::endl;

    std::cout << answer << std::endl;
  }
}