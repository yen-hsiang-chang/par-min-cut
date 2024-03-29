#pragma once

#include <utility>

#include "gbbs/macros.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/random.h"
#include "../utils/utils.h"

#include "cluster.h"

template <class T, class W>
class RCTree {
public:
  RCTree(T n, const parlay::sequence<std::pair<T, T>>& edge_list);
  void build(parlay::random_generator& gen, const parlay::sequence<W>& vertex_value,
             const parlay::sequence<W>& edge_value);
  virtual W f_nullary(Cluster<T, W> *c) = 0;
  virtual W f_unary(Cluster<T, W> *c) = 0;
  virtual W f_binary(Cluster<T, W> *c) = 0;
  virtual W vertex_base() = 0;
  virtual W edge_base() = 0;
  void reevaluate(Cluster<T, W> *c);

protected:
  T n, n_binary, root;
  parlay::sequence<T> child_ptr_binary;
  parlay::sequence<std::pair<T, T>> child_edges_binary, parent_binary;
  parlay::sequence<Cluster<T, W>> edge_clusters;
  parlay::sequence<Cluster<T, W>> vertex_clusters;
  parlay::sequence<Cluster<T, W>> rc_clusters;
};

template <class T, class W>
RCTree<T, W>::RCTree(T n, const parlay::sequence<std::pair<T, T>>& edge_list) : n(n) {
  T n_ternary;
  parlay::sequence<T> edge_ptr(n + 1), edge_ptr_ternary;
  parlay::sequence<std::tuple<T, T, T>> edges(2 * (n - 1)), edges_ternary;

  parlay::parallel_for(0, edges.size(), [&](const T& i) {
    if (i < n - 1)
      edges[i] = std::make_tuple(edge_list[i].first, edge_list[i].second, i);
    else
      edges[i] = std::make_tuple(edge_list[i - (n - 1)].second, edge_list[i - (n - 1)].first, i - (n - 1));
  });

  utils::general_sort(edges, std::less<std::tuple<T, T, T>>());

  parlay::parallel_for(0, edges.size(), [&](const T& i) {
    if (i == 0 || (std::get<0>(edges[i]) != std::get<0>(edges[i - 1]))) {
      edge_ptr[std::get<0>(edges[i])] = i;
    }
  });

  edge_ptr[n] = edges.size();
  
  root = utils::ternary_tree(n, edge_ptr, edges, n_ternary, edge_ptr_ternary, edges_ternary);
  
  utils::binary_tree(root, n_ternary, edge_ptr_ternary, edges_ternary,
                     n_binary, child_ptr_binary, child_edges_binary, parent_binary);
}

template <class T, class W>
void RCTree<T, W>::build(parlay::random_generator& gen, 
                         const parlay::sequence<W>& vertex_value,
                         const parlay::sequence<W>& edge_value) {
  edge_clusters.resize(n_binary - 1);
  parlay::parallel_for(0, n_binary, [&](const T& i) {
    if (i != root) {
      auto [p, id] = parent_binary[i];
      edge_clusters[id] = Cluster<T, W>(id, nullptr, nullptr, nullptr,
                                        nullptr, nullptr, p, i, id, id);
      if (id < n - 1)
        edge_clusters[id].set_val(edge_value[id]);
      else 
        edge_clusters[id].set_val(edge_base());
    }
  });
  
  vertex_clusters.resize(n_binary);
  parlay::parallel_for(0, n_binary, [&](const T& i) {
    vertex_clusters[i] = Cluster<T, W>(i, nullptr, nullptr, nullptr);
    if (i < n)
      vertex_clusters[i].set_val(vertex_value[i]);
    else
      vertex_clusters[i].set_val(vertex_base());
  });

  auto aux_ptr = parlay::sequence<Cluster<T, W>*>::from_function(
    n_binary - 1, [&](const T& i) {return &edge_clusters[i];});
  
  rc_clusters.resize(n_binary);

  auto deg = parlay::sequence<T>::from_function(
    n_binary, [&](const T& i) {return child_ptr_binary[i + 1] - child_ptr_binary[i] + (i != root);});

  parlay::sequence<T> active[2];
  for (T i = 0; i < n_binary; ++i)
    if (i != root && deg[i] > 0 && deg[i] <= 2)
      active[deg[i] - 1].emplace_back(i);

  while(active[0].size() + active[1].size()) {
    gen();
    auto rake = [&](const T& v) {
      Cluster<T, W> *binary_top;
      Cluster<T, W> *unary[2] = {nullptr, nullptr};
      T usz = 0, b0, be0;
      for (T i = child_ptr_binary[v]; i < child_ptr_binary[v + 1]; ++i)
      {
        auto [to, id] = child_edges_binary[i];
        auto &c = aux_ptr[id];
        c -> set_parent_cluster(&rc_clusters[v]);
        unary[usz++] = c;
      }
      {
        auto [to, id] = parent_binary[v];
        auto &c = aux_ptr[id];
        c -> set_parent_cluster(&rc_clusters[v]);
        binary_top = c;
        b0 = (c -> get_top_boundary()) ^ v ^ (c -> get_bottom_boundary());
        be0 = (c -> get_top_boundary_edge()) ^ id ^ (c -> get_bottom_boundary_edge());
        aux_ptr[be0] = &rc_clusters[v];
      }
      rc_clusters[v] = Cluster<T, W>(v, unary[0], unary[1], &vertex_clusters[v],
                                     binary_top, b0, be0);
      vertex_clusters[v].set_parent_cluster(&rc_clusters[v]);
      deg[rc_clusters[v].get_top_boundary()]--;
      deg[v]--;
      rc_clusters[v].set_val(f_unary(&rc_clusters[v]));
    };

    auto compress = [&](const T& v) {
      Cluster<T, W> *binary[2], *unary[2] = {nullptr, nullptr};
      T boundary[2], boundary_edge[2];
      T usz = 0;
      for (T i = child_ptr_binary[v]; i < child_ptr_binary[v + 1]; ++i)
      {
        auto [to, id] = child_edges_binary[i];
        auto &c = aux_ptr[id];
        if (c -> is_binary())
        {
          binary[1] = c;
          boundary[1] = (c -> get_top_boundary()) ^ v ^ (c -> get_bottom_boundary());
          boundary_edge[1] = (c -> get_top_boundary_edge()) ^ id ^ (c -> get_bottom_boundary_edge());
        }
      }
      {
        auto [to, id] = parent_binary[v];
        auto &c = aux_ptr[id];
        binary[0] = c;
        boundary[0] = (c -> get_top_boundary()) ^ v ^ (c -> get_bottom_boundary());
        boundary_edge[0] = (c -> get_top_boundary_edge()) ^ id ^ (c -> get_bottom_boundary_edge());
      }
      if (deg[boundary[1]] <= 1)  return;
      std::bernoulli_distribution dist(0.5);
      auto gen_child = gen[boundary[1]];
      bool color_child = dist(gen_child);
      if (color_child)  return;
      auto gen_v = gen[v];
      bool color_v = dist(gen_v);
      if (!color_v)  return;
      for (T i = child_ptr_binary[v]; i < child_ptr_binary[v + 1]; ++i)
      {
        auto [to, id] = child_edges_binary[i];
        auto &c = aux_ptr[id];
        c -> set_parent_cluster(&rc_clusters[v]);
        if (!(c -> is_binary()))
          unary[usz++] =c;
      }
      {
        auto [to, id] = parent_binary[v];
        auto &c = aux_ptr[id];
        c -> set_parent_cluster(&rc_clusters[v]);
      }
      aux_ptr[boundary_edge[0]] = &rc_clusters[v];
      aux_ptr[boundary_edge[1]] = &rc_clusters[v];
      rc_clusters[v] = Cluster<T, W>(v, unary[0], unary[1], &vertex_clusters[v],
                                     binary[0], binary[1], boundary[0], boundary[1],
                                     boundary_edge[0], boundary_edge[1]);
      vertex_clusters[v].set_parent_cluster(&rc_clusters[v]);
      deg[v] -= 2;
      rc_clusters[v].set_val(f_binary(&rc_clusters[v]));
    };

    for (auto v : active[1])
        compress(v);

    for (auto v : active[0])
        rake(v);

    active[0].clear();
    active[1].clear();
    for (T i = 0; i < n_binary; ++i)
      if (i != root && deg[i] > 0 && deg[i] <= 2)
        active[deg[i] - 1].emplace_back(i);
  }

  auto finalize = [&](const T& v) {
    Cluster<T, W> *unary[2] = {nullptr, nullptr};
    T usz = 0;
    for (T i = child_ptr_binary[v]; i < child_ptr_binary[v + 1]; ++i) {
      auto [to, id] = child_edges_binary[i];
      auto &c = aux_ptr[id];
      c -> set_parent_cluster(&rc_clusters[v]);
      unary[usz++] = c;
    }
    rc_clusters[v] = Cluster(v, unary[0], unary[1], &vertex_clusters[v]);
    vertex_clusters[v].set_parent_cluster(&rc_clusters[v]);
    deg[v] = 0;
    rc_clusters[v].set_val(f_nullary(&rc_clusters[v]));
  };

  finalize(root);
  parlay::parallel_for(0, rc_clusters.size(), [&](const T& i) {
    rc_clusters[i].walk();
  });
  parlay::parallel_for(0, vertex_clusters.size(), [&](const T& i) {
    vertex_clusters[i].walk_from_parent_cluster();
  });
  parlay::parallel_for(0, edge_clusters.size(), [&](const T& i) {
    edge_clusters[i].walk_from_parent_cluster();
  });
  // rc_clusters[root].print();
}

template <class T, class W>
void RCTree<T, W>::reevaluate(Cluster<T, W> *c) {
  while (c = c -> get_parent_cluster()) {
    switch(c -> get_cluster_type()) {
      case 0: c -> set_val(f_nullary(c)); break;
      case 1: c -> set_val(f_unary(c)); break;
      case 2: c -> set_val(f_binary(c)); break;
      default: assert(false);
    }
  }
}