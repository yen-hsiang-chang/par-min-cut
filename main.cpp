#include <iostream>
#include <random>
#include <omp.h>

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/sequence.h"
#include "pbbsbench/benchmarks/comparisonSort/ips4o/sort.h"
#include "gbbs/gbbs/gbbs.h"
#include "gbbs/gbbs/io.h"
#include "gbbs/macros.h"
#include "gbbs/benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"

int main(int argc, char* argv[]) {

  // random number generator
  const size_t seed = 42;
  parlay::random_generator gen(seed);

  // macros
  // using W = double;
  using W = unsigned int;

  // read the graph
  auto G = gbbs::gbbs_io::read_weighted_symmetric_graph<W>(argv[1], false, false);
    // "/global/cfs/cdirs/mp156/yenhsiang/rhg.txt", false, false);
    // "/global/cfs/cdirs/mp156/yenhsiang/test.txt", false, false);
    // "/global/cfs/cdirs/mp156/yenhsiang/largerhg.txt", false, false);
    // "/global/cfs/cdirs/mp156/yenhsiang/superlargerhg.txt", false, false);
    // "/global/cfs/cdirs/mp156/yenhsiang/livejournal.txt", false, false);
  // auto edges = G.edges();
  // auto edges_with_exp = parlay::tabulate(G.m, [&](const size_t& i) {
  //   return std::make_pair(edges[i], double(0));
  // });

  // for(int iter = 0;iter < 1;++iter) {
  //   // Step 1: Generate weighted random edge ordering
  //   auto weighted_ordering_begin = omp_get_wtime();
  //   gen();
  //   parlay::parallel_for(0, G.m, [&](const size_t& i) {
  //     std::exponential_distribution<double> exponential(std::get<2>(edges_with_exp[i].first));
  //     auto local_rand = gen[i];
  //     edges_with_exp[i].second = exponential(local_rand);
  //   });

  //   auto pred = [&](const std::pair<std::tuple<gbbs::uintE, gbbs::uintE, W>, double>& l,
  //                   const std::pair<std::tuple<gbbs::uintE, gbbs::uintE, W>, double>& r) {
  //                     return l.second < r.second;
  //                   };
  //   compSort(edges_with_exp, pred);
  //   auto weighted_ordering_end = omp_get_wtime();
  //   std::cout << "Weighted Random Edge Ordering Time: "
  //             << double(weighted_ordering_end - weighted_ordering_begin) << " seconds" << std::endl;
  
    // Step 2: Compute MST based on weighted random edge ordering
    auto mst_begin = omp_get_wtime();
    gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G);
    auto mst_end = omp_get_wtime();

    std::cout << "MST Time: "
              << double(mst_end - mst_begin) << " seconds" << std::endl;
  // }
}