#include <iostream>
#include <random>
#include <omp.h>

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "gbbs/gbbs/gbbs.h"
#include "gbbs/gbbs/io.h"
#include "gbbs/macros.h"
#include "gbbs/benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "rctree/rctree.h"
#include "utils/utils.h"

int main(int argc, char* argv[]) {

  // random number generator
  const size_t seed = 42;
  parlay::random_generator gen(seed);

  // macros
  using W = double;
  // using W = unsigned int;

  // read the graph, assume {u, v} only has (u, v) or (v, u)
  auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<W>(argv[1], false, false);

  // (v, Exp(w))
  auto E_MST_input = gbbs::new_array_no_init<std::tuple<gbbs::uintE, gbbs::uintE>>(G.m);
  parlay::parallel_for(0, G.m, [&](const size_t& i) {
    E_MST_input[i] = std::make_tuple(std::get<0>(G.e0[i]), 0);
  });

  auto G_MST_input = gbbs::symmetric_graph<gbbs::symmetric_vertex, gbbs::uintE>(
      G.v_data, G.n, G.m,
      [](){},
      E_MST_input);
  
  auto E_exponential = parlay::sequence<std::pair<double, gbbs::uintE>>(G.m);

  for(int iter = 0;iter < 1;++iter) {
    // Step 1: Generate weighted random edge ordering
    auto weighted_ordering_begin = omp_get_wtime();
    gen();
    parlay::parallel_for(0, G.m, [&](const size_t& i) {
      std::exponential_distribution<double> exponential(std::get<1>(G.e0[i]));
      auto local_rand = gen[i];
      E_exponential[i] = std::make_pair(exponential(local_rand), gbbs::uintE(i));
    });

    general_sort(E_exponential, std::less<std::pair<double, gbbs::uintE>>());

    parlay::parallel_for(0, G.m, [&](const size_t& i) {
      std::get<1>(E_MST_input[E_exponential[i].second]) = i; 
    });
      
    auto weighted_ordering_end = omp_get_wtime();
    std::cout << "Weighted Random Edge Ordering Time: "
              << double(weighted_ordering_end - weighted_ordering_begin) << " seconds" << std::endl;
  
    // Step 2: Compute MST based on weighted random edge ordering
    auto mst_begin = omp_get_wtime();
    auto E_MST = gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G_MST_input);
    auto mst_end = omp_get_wtime();

    std::cout << "MST Time: "
              << double(mst_end - mst_begin) << " seconds" << std::endl;
    
    if (E_MST.size() + 1 != G.n)
    {
      std::cout << "min cut = 0." << std::endl; 
      break;
    }

    auto MST_edge_list = parlay::sequence<std::pair<gbbs::uintE, gbbs::uintE>>::from_function(
      E_MST.size(), [&](size_t i) {
        auto [u, v, _] = E_MST[i];
        return std::make_pair(u, v);
      });
    
    auto vertex_weight = parlay::sequence<W>(G.n, 1);
    auto edge_weight = parlay::sequence<W>(G.n - 1, 2);

    auto rctree = RCTree<W>(G.n, MST_edge_list, vertex_weight, edge_weight);

  }

  gbbs::free_array(E_MST_input, G.m);
}