#pragma once

#include <utility>

#include "gbbs/macros.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "../utils/utils.h"

#include "cluster.h"

template<class T>
class RCTree {
public:
  RCTree(gbbs::uintE n, parlay::sequence<std::pair<gbbs::uintE, gbbs::uintE>> edge_list, 
               parlay::sequence<T> vertex_weight, parlay::sequence<T> edge_weight);

private:
  gbbs::uintE n;
  size_t m;
  parlay::sequence<size_t> edge_ptr;
  parlay::sequence<Binary_Cluster> edge_clusters;
  parlay::sequence<Nullary_Cluster> vertex_clusters;
  parlay::sequence<std::tuple<gbbs::uintE, gbbs::uintE, size_t>> edges;
};

template<class T>
RCTree<T>::RCTree(gbbs::uintE n, parlay::sequence<std::pair<gbbs::uintE, gbbs::uintE>> edge_list, 
               parlay::sequence<T> vertex_weight, parlay::sequence<T> edge_weight) 
               : n(n), m(2LL * (n - 1)), edge_ptr(n + 1, 0), edges(2 * (n - 1)) {
  
  assert(edge_list.size() == n - 1);
  parlay::parallel_for(0, m, [&](const size_t& i) {
    if (i < n - 1)
      edges[i] = std::make_tuple(edge_list[i].first, edge_list[i].second, i);
    else
      edges[i] = std::make_tuple(edge_list[i - (n - 1)].second, edge_list[i - (n - 1)].first, i - (n - 1));
  });
  general_sort(edges, std::less<std::tuple<gbbs::uintE, gbbs::uintE, size_t>>());
  parlay::parallel_for(0, m, [&](size_t i) {
    if (i == 0 || (std::get<0>(edges[i]) != std::get<0>(edges[i - 1]))) {
      edge_ptr[std::get<0>(edges[i])] = i;
    }
  });
  edge_ptr[n] = m;
}