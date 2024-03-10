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
  RCTree(uint n, parlay::sequence<std::pair<uint, uint>> edge_list, 
         parlay::sequence<T> vertex_weight, parlay::sequence<T> edge_weight,
         parlay::random_generator& gen);

private:
  uint n, m;
  uint n_ternary, m_ternary;
  parlay::sequence<uint> edge_ptr, edge_ptr_ternary;
  parlay::sequence<std::tuple<uint, uint, uint>> edges, edges_ternary;
  parlay::sequence<Cluster> edge_clusters;
  parlay::sequence<Cluster> vertex_clusters;
  parlay::sequence<Cluster> rc_clusters;
};

template<class T>
RCTree<T>::RCTree(uint n, parlay::sequence<std::pair<uint, uint>> edge_list, 
                  parlay::sequence<T> vertex_weight, parlay::sequence<T> edge_weight,
                  parlay::random_generator& gen) 
                  : n(n), m(2 * (n - 1)), edge_ptr(n + 1), edges(2 * (n - 1)) {

  assert(edge_list.size() == n - 1);

  parlay::parallel_for(0, m, [&](const uint& i) {
    if (i < n - 1)
      edges[i] = std::make_tuple(edge_list[i].first, edge_list[i].second, i);
    else
      edges[i] = std::make_tuple(edge_list[i - (n - 1)].second, edge_list[i - (n - 1)].first, i - (n - 1));
  });

  utils::general_sort(edges, std::less<std::tuple<uint, uint, uint>>());

  parlay::parallel_for(0, m, [&](const uint& i) {
    if (i == 0 || (std::get<0>(edges[i]) != std::get<0>(edges[i - 1]))) {
      edge_ptr[std::get<0>(edges[i])] = i;
    }
  });

  edge_ptr[n] = m;

  auto deg = parlay::sequence<uint>::from_function(
    n, [&](uint i) {return edge_ptr[i + 1] - edge_ptr[i];});
  
  parlay::sequence<uint> deg_ternary;
  
  utils::ternary_tree(n, m, edge_ptr, edges, deg,
                      n_ternary, m_ternary, edge_ptr_ternary, edges_ternary, deg_ternary);
  
  edge_clusters = parlay::sequence<Cluster>::from_function(
    n - 1, [&](uint i) {
      return Cluster(i, parlay::sequence<Cluster*>(), nullptr,
              nullptr, nullptr, edge_list[i].first, edge_list[i].second, i, i);
    }
  );
  
  vertex_clusters = parlay::sequence<Cluster>::from_function(
    n, [&](uint i) {
      return Cluster(i, parlay::sequence<Cluster*>(), nullptr);
    }
  );

  auto aux_ptr = parlay::sequence<Cluster*>::from_function(
    n - 1, [&](uint i) {return &edge_clusters[i];});
  
  rc_clusters.resize(n);

  parlay::sequence<uint> active;
  for (uint i = 1;i < n; ++i)
    if (deg[i] == 1 || deg[i] == 2)
      active.emplace_back(i);

  while(!active.empty()) {
    gen();
    auto rake = [&](const uint& v) {
      Cluster *binary_top;
      parlay::sequence<Cluster*> unary;
      uint b0, be0;
      for (uint i = edge_ptr[v]; i < edge_ptr[v + 1]; ++i)
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

    auto compress = [&](const uint& v) {
      Cluster *binary_clusters[2];
      uint boundary[2], boundary_edge[2];
      int bsz = 0;
      for (uint i = edge_ptr[v]; i < edge_ptr[v + 1]; ++i)
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
      bool color_r = boundary[1] == 0 ? false : dist(gen_r);
      if (color_r)  return;
      parlay::sequence<Cluster*> unary;
      for (uint i = edge_ptr[v]; i < edge_ptr[v + 1]; ++i)
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
    for (uint i = 1;i < n; ++i)
      if (deg[i] == 1 || deg[i] == 2)
        active.emplace_back(i);
  }

  auto finalize = [&](const uint& v) {
      parlay::sequence<Cluster*> unary;
      for (uint i = edge_ptr[v]; i < edge_ptr[v + 1]; ++i)
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