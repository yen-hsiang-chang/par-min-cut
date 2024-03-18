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
  using W = double;

  // random number generator
  const size_t seed = 42;
  parlay::random_generator gen(seed);

  // read the graph, assume {u, v} only has (u, v) or (v, u)
  auto G_dir = gbbs::gbbs_io::read_weighted_symmetric_graph<W>(argv[1], false, false);
  auto G_edge_list = gbbs::to_edge_array<W>(G_dir);
  auto G_undir = gbbs::gbbs_io::edge_list_to_symmetric_graph<W>(G_edge_list);

  // (v, Exp(w))
  auto E_MST_input = parlay::sequence<std::tuple<uint, double>>(G_dir.m);
  auto V_MST_input = parlay::sequence<gbbs::vertex_data>(G_dir.n);

  auto G_MST_input = gbbs::symmetric_graph<gbbs::symmetric_vertex, double>(
    V_MST_input.data(), G_dir.n, G_dir.m, [](){}, E_MST_input.data());
  
  auto contraction_vertex_weight = parlay::sequence<Contraction_Type<W>>(G_undir.n);
  auto G_weighted_degree = parlay::sequence<W>(G_undir.n);

  parlay::parallel_for(0, G_undir.n, [&](const uint& i) {
    auto get_weight = [&](const size_t& j) {
      return std::get<1>(G_undir.e0[G_undir.v_data[i].offset + j]); 
    };
    auto weights = parlay::delayed_seq<W>(G_undir.v_data[i].degree, get_weight);
    G_weighted_degree[i] = parlay::reduce(weights);
    contraction_vertex_weight[i] = Contraction_Type<W>(G_weighted_degree[i], 0, 0, 0, true);
  });
  auto contraction_edge_weight = parlay::sequence<Contraction_Type<double>>::from_function(
    G_undir.n - 1, [&](const uint& i) {
      return Contraction_Type<W>(0, 0, i, i, false);
    });
  W G_min_weighted_degree = parlay::reduce(G_weighted_degree, parlay::minimum<W>());

  for(int iter = 0; iter < 100; ++iter) {
    // Step 1: Generate weighted random edge ordering
    auto weighted_ordering_time = -omp_get_wtime();
    gen();
    parlay::parallel_for(0, G_dir.m, [&](const uint& i) {
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
    // Theses will be overwritten in MST algo
    G_MST_input.n = G_dir.n;
    G_MST_input.m = G_dir.m;
    parlay::parallel_for(0, G_dir.n, [&](const uint& i) {
      V_MST_input[i] =  G_dir.v_data[i];
    });
    auto E_MST = gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G_MST_input);
    mst_time += omp_get_wtime();

    std::cout << "MST Time: "
              << mst_time << " seconds" << std::endl;
    
    if (E_MST.size() + 1 != G_dir.n)
    {
      std::cout << "min cut = 0." << std::endl; 
      break;
    }

    auto contraction_rctree_time = -omp_get_wtime();
    utils::general_sort(E_MST, [&](const std::tuple<uint, uint, double>& lhs, 
                                   const std::tuple<uint, uint, double>& rhs) {
                                    return std::get<2>(lhs) < std::get<2>(rhs);
                                   });

    auto MST_edge_list = parlay::sequence<std::pair<uint, uint>>::from_function(
      E_MST.size(), [&](const uint& i) {
        auto [u, v, _] = E_MST[i];
        return std::make_pair(u, v);
      });

    auto rctree = Contraction_RCTree<uint, W>(G_undir.n, MST_edge_list);
    rctree.build(gen, contraction_vertex_weight, contraction_edge_weight);
    contraction_rctree_time += omp_get_wtime();

    std::cout << "Contraction RCTree Time: "
              << contraction_rctree_time << " seconds" << std::endl;
    auto contraction_mixop = parlay::sequence<Contraction_MixOp<W>>((G_edge_list.size() + MST_edge_list.size()) * 2);
    parlay::parallel_for(0, G_edge_list.E.size(), [&](const size_t& i) {
      auto [u, v, w] = G_edge_list.E[i];
      contraction_mixop[i * 2] = Contraction_MixOp<W>(rctree.query_max_path(u, v) * 3, 0, u, w);
      contraction_mixop[i * 2 + 1] = Contraction_MixOp<W>(rctree.query_max_path(u, v) * 3, 0, v, w);
    });
    parlay::parallel_for(0, MST_edge_list.size(), [&](const size_t& i) {
      auto [u, v] = MST_edge_list[i];
      contraction_mixop[G_edge_list.size() * 2 + i * 2] = Contraction_MixOp<W>(i * 3 + 1, 1, u, 0);
      contraction_mixop[G_edge_list.size() * 2 + i * 2 + 1] = Contraction_MixOp<W>(i * 3 + 2, 2, u, 0);
    });
    utils::general_sort(contraction_mixop, [&](const Contraction_MixOp<W>& lhs, 
                                               const Contraction_MixOp<W>& rhs) {
                                                return lhs.timestamp < rhs.timestamp;
                                               });
    W answer = G_min_weighted_degree;
    for (auto q : contraction_mixop) {
      if (q.type == 0) {
        rctree.subtract_vertex_weight(q.u, q.w);
      } else if(q.type == 1) {
        rctree.join_edge(q.timestamp / 3);
      } else {
        if (q.timestamp != 3 * MST_edge_list.size() - 1)
          answer = std::min(answer, rctree.query_vertex_sum(q.u));
      }
    }
    std::cout << answer << std::endl;
  }
}