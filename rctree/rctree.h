#pragma once

#include <utility>

#include "gbbs/macros.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/random.h"
#include "../utils/utils.h"

#include "cluster.h"

template<class T>
class RCTree {
public:
  RCTree(gbbs::uintE n, parlay::sequence<std::pair<gbbs::uintE, gbbs::uintE>> edge_list, 
         parlay::sequence<T> vertex_weight, parlay::sequence<T> edge_weight,
         parlay::random_generator& gen);

private:
  gbbs::uintE n;
  size_t m;
  parlay::sequence<size_t> edge_ptr;
  parlay::sequence<std::tuple<gbbs::uintE, gbbs::uintE, size_t>> edges;
  parlay::sequence<Cluster> edge_clusters;
  parlay::sequence<Cluster> vertex_clusters;
  parlay::sequence<Cluster> rc_clusters;
};

template<class T>
RCTree<T>::RCTree(gbbs::uintE n, parlay::sequence<std::pair<gbbs::uintE, gbbs::uintE>> edge_list, 
                  parlay::sequence<T> vertex_weight, parlay::sequence<T> edge_weight,
                  parlay::random_generator& gen) 
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

  auto deg = parlay::sequence<size_t>::from_function(
    n, [&](gbbs::uintE i) {return edge_ptr[i + 1] - edge_ptr[i];});
  
  edge_clusters = parlay::sequence<Cluster>::from_function(
    n - 1, [&](gbbs::uintE i) {
      return Cluster(i, parlay::sequence<Cluster*>(), nullptr,
              nullptr, nullptr, edge_list[i].first, edge_list[i].second, i, i);
    }
  );
  
  vertex_clusters = parlay::sequence<Cluster>::from_function(
    n, [&](gbbs::uintE i) {
      return Cluster(i, parlay::sequence<Cluster*>(), nullptr);
    }
  );

  auto aux_ptr = parlay::sequence<Cluster*>::from_function(
    n - 1, [&](gbbs::uintE i) {return &edge_clusters[i];});
  
  rc_clusters.resize(n);

  parlay::sequence<gbbs::uintE> active;
  for (gbbs::uintE i = 1;i < n; ++i)
    if (deg[i] == 1 || deg[i] == 2)
      active.emplace_back(i);

  while(!active.empty()) {
    gen();
    auto rake = [&](const gbbs::uintE& v) {
      Cluster *binary_top;
      parlay::sequence<Cluster*> unary;
      gbbs::uintE b0, be0;
      for (gbbs::uintE i = edge_ptr[v]; i < edge_ptr[v + 1]; ++i)
      {
        auto [from, to, id] = edges[i];
        auto &c = aux_ptr[id];
        c -> parent = &rc_clusters[v];
        if (c -> is_binary())
        {
          binary_top = c;
          b0 = (c -> boundary[0]) ^ v ^ (c -> boundary[1]);
          be0 = (c -> boundary_edge[0]) ^ id ^ (c -> boundary_edge[1]);
          aux_ptr[be0] = &rc_clusters[v];
        }
        else
          unary.emplace_back(c);
      }
      rc_clusters[v] = Cluster(v, unary, &vertex_clusters[v],
                               binary_top, b0, be0);
      vertex_clusters[v].parent = &rc_clusters[v];
      deg[rc_clusters[v].boundary[0]]--;
      deg[v]--;
    };

    auto compress = [&](const gbbs::uintE& v) {
      Cluster *binary_clusters[2];
      gbbs::uintE boundary[2], boundary_edge[2];
      int bsz = 0;
      for (gbbs::uintE i = edge_ptr[v]; i < edge_ptr[v + 1]; ++i)
      {
        auto [from, to, id] = edges[i];
        auto &c = aux_ptr[id];
        if (c -> is_binary())
        {
          binary_clusters[bsz] = c;
          boundary[bsz] = (c -> boundary[0]) ^ v ^ (c -> boundary[1]);
          boundary_edge[bsz] = (c -> boundary_edge[0]) ^ id ^ (c -> boundary_edge[1]);
          bsz++;
        }
      }
      if (boundary[0] != 0 && deg[boundary[0]] <= 1)  return;
      if (boundary[1] != 0 && deg[boundary[1]] <= 1)  return;
      std::bernoulli_distribution dist(0.5);
      auto gen_l = gen[boundary[0]];
      bool color_l = boundary[0] == 0 ? false : dist(gen_l);
      if (color_l)  return;
      auto gen_m = gen[v];
      bool color_m = dist(gen_m);
      if (!color_m)  return;
      auto gen_r = gen[boundary[1]];
      bool color_r = boundary[0] == 0 ? false : dist(gen_r);
      if (color_r)  return;
      parlay::sequence<Cluster*> unary;
      for (gbbs::uintE i = edge_ptr[v]; i < edge_ptr[v + 1]; ++i)
      {
        auto [from, to, id] = edges[i];
        auto &c = aux_ptr[id];
        c -> parent = &rc_clusters[v];
        if (!(c -> is_binary()))
          unary.emplace_back(c);
      }
      aux_ptr[boundary_edge[0]] = &rc_clusters[v];
      aux_ptr[boundary_edge[1]] = &rc_clusters[v];
      rc_clusters[v] = Cluster(v, unary, &vertex_clusters[v], binary_clusters[0], binary_clusters[1],
                               boundary[0], boundary[1], boundary_edge[0], boundary_edge[1]);
      vertex_clusters[v].parent = &rc_clusters[v];
      deg[v] -= 2;
    };

    for (auto v : active) {
      if (deg[v] == 2)
        compress(v);
    }

    for (auto v : active) {
      if (deg[v] == 1)
        rake(v);
    }

    active.clear();
    for (gbbs::uintE i = 1;i < n; ++i)
      if (deg[i] == 1 || deg[i] == 2)
        active.emplace_back(i);
  }

  auto finalize = [&](const gbbs::uintE& v) {
      parlay::sequence<Cluster*> unary;
      for (gbbs::uintE i = edge_ptr[v]; i < edge_ptr[v + 1]; ++i)
      {
        auto [from, to, id] = edges[i];
        auto &c = aux_ptr[id];
        c -> parent = &rc_clusters[v];
        unary.emplace_back(c);
      }
      rc_clusters[v] = Cluster(v, unary, &vertex_clusters[v]);
      vertex_clusters[v].parent = &rc_clusters[v];
      deg[v] = 0;
    };

  finalize(0);
  rc_clusters[0].print();
}