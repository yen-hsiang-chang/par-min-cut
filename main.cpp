#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <omp.h>
#include <unistd.h>

#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/random.h"
#include "gbbs/gbbs/gbbs.h"
#include "gbbs/gbbs/graph_io.h"
#include "gbbs/gbbs/io.h"
#include "gbbs/macros.h"
#include "gbbs/benchmarks/MinimumSpanningForest/Boruvka/MinimumSpanningForest.h"
#include "rctree/contraction/contraction.h"
#include "utils/utils.h"

int main(int argc, char* argv[]) {
  // macros
  using W = double;
  int opt = -1;
  std::string src_graph = "";
  size_t num_iter = 1;
  size_t seed = 42;
  enum Paradigm {parrec, seqrec_uf, seqrec_rc};
  Paradigm paradigm = parrec;
  
  auto usage = []() {
    std::cerr << "\nUsage: ./build/min-cut [options]\n\n"
              << "Options:\n"
              << "    -g <Source Graph File Path>                 Path of file with source graph in edge list format.\n"
              << "    -k <Number of rounds of edge contractions>  Default to 1, check the referenced paper to see how many rounds are needed to achieve the desired approximation w.h.p.\n"
              << "    -s <Seed>                                   Default to 42, the seed used for randomized algorithms.\n"
              << "    -p <Paradigm>                               Choose among parrec (default), seqrec-uf, and seqrec-rc. Check the report for more information.\n"
              << "    -h                                          Help.\n";
  };
  while((opt = getopt(argc, argv, "g:k:s:p:h")) != -1) {
    switch (opt) {
      case 'g':
        src_graph = optarg;
        break;
      case 'k':
        num_iter = std::stoll(optarg);
        break;
      case 's':
        seed = std::stoll(optarg);
        break;
      case 'p':
        if (std::string(optarg) == "parrec") {
          paradigm = parrec;
        } else if (std::string(optarg) == "seqrec-uf") {
          paradigm = seqrec_uf;
        } else if (std::string(optarg) == "seqrec-rc") {
          paradigm = seqrec_rc;
        } else {
          std::cerr << "Paradigm not supported!\n";
          usage();
          exit(0);
        }
        break;
      case 'h':
        usage();
        exit(0);
        break;
      default:
        std::cerr << "Unrecognized option!\n";
        usage();
        exit(0);
    }
  }

  if (src_graph == "") {
    std::cerr << "Source graph is not provided!\n";
    usage();
    exit(0);
  }

  // random number generator
  parlay::random_generator gen(seed);

  // read the graph
  auto G_edge_list = gbbs::gbbs_io::read_weighted_edge_list<W>(src_graph.c_str());
  auto G_undir = gbbs::gbbs_io::edge_list_to_symmetric_graph<W>(G_edge_list);

  // Feed exponential distribution to MST solver with (v, Exp(w))
  auto E_MST_input = parlay::sequence<std::tuple<uintV, double>>(G_undir.m);
  auto V_MST_input = parlay::sequence<gbbs::vertex_data>(G_undir.n);
  auto G_MST_input = gbbs::symmetric_graph<gbbs::symmetric_vertex, double>(
    V_MST_input.data(), G_undir.n, G_undir.m, [](){}, E_MST_input.data());
  
  // Construct vertex weight for Contraction RCTree
  auto contraction_vertex_weight = parlay::sequence<Contraction_Type<W>>(G_undir.n);
  auto G_weighted_degree = parlay::sequence<W>(G_undir.n);
  parlay::parallel_for(0, G_undir.n, [&](uintV i) {
    auto get_weight = [&](uintV j) {
      return std::get<1>(G_undir.e0[G_undir.v_data[i].offset + j]); 
    };
    auto weights = parlay::delayed_seq<W>(G_undir.v_data[i].degree, get_weight);
    G_weighted_degree[i] = parlay::reduce(weights);
    contraction_vertex_weight[i] = std::make_pair(Contraction_CC<W>(G_weighted_degree[i], 0, true), Contraction_P(0, 0));
  });

  // Construct edge weight for Contraction RCTree
  auto contraction_edge_weight = parlay::sequence<Contraction_Type<W>>::from_function(G_undir.n - 1, [&](uintV i) {
    return std::make_pair(Contraction_CC<W>(0, 0, false), Contraction_P(i, i)); 
  });

  // Base case where each vertex is a partition
  W G_min_weighted_degree = parlay::reduce(G_weighted_degree, parlay::minimum<W>()), global_mc = std::numeric_limits<W>::max();

  auto total_time = -omp_get_wtime();
  
  // TODO: determine the number of iterations from n
  for(size_t iter = 0; iter < num_iter; ++iter) {
    std::cout << "\nIteration: " << iter << '\n';
    // Step 1: Generate weighted random edge ordering via exponential distribution
    auto weighted_ordering_time = -omp_get_wtime();
    gen(); // Advance one step for the rng
    parlay::parallel_for(0, G_undir.m, [&](uintE i) {
      auto [v, w] = G_undir.e0[i];
      std::exponential_distribution<double> exponential(w);
      auto local_rand = gen[i];
      E_MST_input[i] = std::make_tuple(v, exponential(local_rand));
    });
    weighted_ordering_time += omp_get_wtime();

    std::cout << "    Weighted Random Edge Ordering Time: "
              << weighted_ordering_time << " seconds" << "\n";
  
    // Step 2: Compute MST based on weighted random edge ordering
    auto mst_time = -omp_get_wtime();
    // These will be overwritten in MST algo, reset every iteration
    G_MST_input.n = G_undir.n;
    G_MST_input.m = G_undir.m;
    parlay::parallel_for(0, G_undir.n, [&](uintV i) {
      V_MST_input[i] =  G_undir.v_data[i];
    });
    auto E_MST = gbbs::MinimumSpanningForest_boruvka::MinimumSpanningForest(G_MST_input);
    mst_time += omp_get_wtime();

    std::cout << "    MST Time: "
              << mst_time << " seconds" << "\n";
    
    // If it is a forest, the minimum cut is just 0.
    if (E_MST.size() + 1 != G_undir.n) {
      std::cout << "min cut = 0." << "\n"; 
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

    std::cout << "    Contraction RCTree Time: "
              << contraction_rctree_time << " seconds" << "\n";

    // Step 4: Batch operations
    if (paradigm != seqrec_uf) {
      auto rec_time = -omp_get_wtime();
      auto contraction_mixop = parlay::sequence<Contraction_MixOp<W>>((G_edge_list.size() + MST_edge_list.size() - 1) * 2);
      // Subtract vertex weight
      parlay::parallel_for(0, G_edge_list.size(), [&](uintE i) {
        auto [u, v, w] = G_edge_list[i];
        contraction_mixop[i * 2] = Contraction_MixOp<W>(contraction_rctree.query_max_path(u, v) * 4 + 1, 0, u, w);
        contraction_mixop[i * 2 + 1] = Contraction_MixOp<W>(contraction_rctree.query_max_path(u, v) * 4 + 2, 0, v, w);
      });
      // Join edges and Query component induced by joined edges
      parlay::parallel_for(0, MST_edge_list.size() - 1, [&](uintE i) {
        auto [u, v] = MST_edge_list[i];
        contraction_mixop[G_edge_list.size() * 2 + i * 2] = Contraction_MixOp<W>(i * 4 + 3, 1, i, 0);
        contraction_mixop[G_edge_list.size() * 2 + i * 2 + 1] = Contraction_MixOp<W>(i * 4 + 4, 2, u, 0);
      });
      W answer = paradigm == parrec ? contraction_rctree.batch_operations(contraction_mixop) : 
                                      contraction_rctree.batch_operations_sequential(contraction_mixop);
      rec_time += omp_get_wtime();

      if (paradigm == parrec) {
        std::cout << "    ParREC Time: "
                  << rec_time << " seconds" << "\n";
      } else {
        std::cout << "    SeqREC-RC Time: "
                  << rec_time << " seconds" << "\n";
      }

      std::cout << "    Approximate Minimum Cut: " << answer << "\n";
      global_mc = std::min(global_mc, answer);
    } else {
      std::vector<uintV> p(G_undir.n);
      std::vector<W> sz(G_undir.n);
      for (uintV i = 0;i < G_undir.n; ++i)
        sz[i] = G_weighted_degree[i], p[i] = i;
      std::function<uintV(uintV)> F = [&](uintV x) {
        return x == p[x] ? x : p[x] = F(p[x]);
      };
      std::function<void(uintV, uintV)> U = [&](uintV x, uintV y) {
        x = F(x), y = F(y);
        if(x != y) {
          p[x] = y, sz[y] += sz[x];
        }
      };
      std::function<void(uintV, W)> D = [&](uintV x, W w) {
        sz[F(x)] -= w;
      };
      std::function<W(uintV)> Q = [&](uintV x) {
        return sz[F(x)];
      };
      auto seqrec_uf_time = -omp_get_wtime();
      auto contraction_mixop = parlay::sequence<Contraction_MixOp<W>>((G_edge_list.size() + MST_edge_list.size() - 1) * 2);
      // Subtract vertex weight
      parlay::parallel_for(0, G_edge_list.size(), [&](uintE i) {
        auto [u, v, w] = G_edge_list[i];
        contraction_mixop[i * 2] = Contraction_MixOp<W>(contraction_rctree.query_max_path(u, v) * 4 + 1, 0, u, w);
        contraction_mixop[i * 2 + 1] = Contraction_MixOp<W>(contraction_rctree.query_max_path(u, v) * 4 + 2, 0, v, w);
      });
      // Join edges and Query component induced by joined edges
      parlay::parallel_for(0, MST_edge_list.size() - 1, [&](uintE i) {
        auto [u, v] = MST_edge_list[i];
        contraction_mixop[G_edge_list.size() * 2 + i * 2] = Contraction_MixOp<W>(i * 4 + 3, 1, i, 0);
        contraction_mixop[G_edge_list.size() * 2 + i * 2 + 1] = Contraction_MixOp<W>(i * 4 + 4, 2, u, 0);
      });
      utils::general_sort(contraction_mixop, [&](const Contraction_MixOp<W>& lhs, const Contraction_MixOp<W>& rhs) {
        return lhs.timestamp < rhs.timestamp;
      });
      W answer = std::numeric_limits<W>::max();
      for (auto q : contraction_mixop) {
        if (q.type == 0) {
          D(q.u, q.w);
        } else if(q.type == 1) {
          U(MST_edge_list[q.u].first, MST_edge_list[q.u].second);
        } else {
          answer = std::min(answer, W(Q(q.u)));
        }
      }
      seqrec_uf_time += omp_get_wtime();

      std::cout << "    SeqREC-UF Time: "
                << seqrec_uf_time << " seconds" << "\n";
      
      std::cout << "    Approximate Minimum Cut: " << answer << "\n";
      global_mc = std::min(global_mc, answer);
    }
  }

  total_time += omp_get_wtime();

  std::cout << "\nTotal Time: "
            << total_time << " seconds" << "\n"
            << "Overall Appoximate Minimum Cut: " << std::min(G_min_weighted_degree, global_mc) << "\n";
}