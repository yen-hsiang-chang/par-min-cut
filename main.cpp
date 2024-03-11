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
#include "rctree/rctree.h"
#include "utils/utils.h"

int main(int argc, char* argv[]) {

  using namespace utils;

  // random number generator
  const size_t seed = 42;
  parlay::random_generator gen(seed);

  // macros
  using W = double;

  // read the graph, assume {u, v} only has (u, v) or (v, u)
  auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<W>(argv[1], false, false);

  // (v, Exp(w))
  auto E_MST_input = gbbs::new_array_no_init<std::tuple<uint, double>>(G.m);

  auto G_MST_input = gbbs::symmetric_graph<gbbs::symmetric_vertex, double>(
      G.v_data, G.n, G.m, [](){}, E_MST_input);

  for(int iter = 0;iter < 1;++iter) {
    // Step 1: Generate weighted random edge ordering
    auto weighted_ordering_begin = omp_get_wtime();
    gen();
    parlay::parallel_for(0, G.m, [&](const size_t& i) {
      auto [v, w] = G.e0[i];
      std::exponential_distribution<double> exponential(w);
      auto local_rand = gen[i];
      E_MST_input[i] = std::make_tuple(v, exponential(local_rand));
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

    auto MST_edge_list = parlay::sequence<std::pair<uint, uint>>::from_function(
      E_MST.size(), [&](size_t i) {
        auto [u, v, _] = E_MST[i];
        return std::make_pair(u, v);
      });
    
    // auto vertex_weight = parlay::sequence<W>(G.n, 1);
    // auto edge_weight = parlay::sequence<W>(G.n - 1, 2);

    auto rctree = RCTree<uint>(G.n, MST_edge_list);

  }

  gbbs::free_array(E_MST_input, G.m);
}