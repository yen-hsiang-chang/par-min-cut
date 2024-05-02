#pragma once

#include <utility>

#include "gbbs/macros.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/random.h"
#include "../utils/utils.h"

#include "cluster.h"

template <class W>
class RCTree {
public:
  // Constructor where build needs to be called in the constructor of inherited class
  // since it needs the pure virtual functions.
  RCTree(uintV n, const parlay::sequence<std::pair<uintV, uintV>>& edge_list);
  void build(parlay::random_generator& gen, const parlay::sequence<W>& vertex_value,
             const parlay::sequence<W>& edge_value);
  
  // Pure Virtual functions that need to be defined in inherited class for maintaining values in nullary/unary/binary clusters
  virtual W f_nullary(Cluster<W> *c) = 0;
  virtual W f_unary(Cluster<W> *c) = 0;
  virtual W f_binary(Cluster<W> *c) = 0;

  // Pure Virtual functions that need to be defined in inherited class for added vertex/edge
  virtual W vertex_base() = 0;
  virtual W edge_base() = 0;

protected:
  // Number of vertices in the original tree
  uintV n;
  // Data structure for the rooted binary tree
  uintV n_binary, root;
  parlay::sequence<uintV> child_ptr_binary;
  parlay::sequence<std::pair<uintV, uintV>> child_edges_binary, parent_binary;
  // Clusters in RCTree
  parlay::sequence<Cluster<W>> edge_clusters;
  parlay::sequence<Cluster<W>> vertex_clusters;
  parlay::sequence<Cluster<W>> rc_clusters;

  // Walking up the tree from leaf cluster to reevaluate value
  void reevaluate(Cluster<W> *c);
};

// Constructor for RCTree where it first converts the tree to a rooted binary tree
template <class W>
RCTree<W>::RCTree(uintV n, const parlay::sequence<std::pair<uintV, uintV>>& edge_list) : n(n) {
  uintV n_ternary;
  parlay::sequence<uintV> edge_ptr(n + 1);
  parlay::sequence<uint64_t> edge_ptr_ternary;
  parlay::sequence<std::tuple<uintV, uintV, uintV>> edges(2 * (n - 1)), edges_ternary;

  parlay::parallel_for(0, edges.size(), [&](uintV i) {
    if (i < n - 1)
      edges[i] = std::make_tuple(edge_list[i].first, edge_list[i].second, i);
    else
      edges[i] = std::make_tuple(edge_list[i - (n - 1)].second, edge_list[i - (n - 1)].first, i - (n - 1));
  });

  utils::general_sort(edges, std::less<std::tuple<uintV, uintV, uintV>>());

  parlay::parallel_for(0, edges.size(), [&](uintV i) {
    if (i == 0 || (std::get<0>(edges[i]) != std::get<0>(edges[i - 1]))) {
      edge_ptr[std::get<0>(edges[i])] = i;
    }
  });

  edge_ptr[n] = edges.size();
  
  root = utils::ternary_tree(n, edge_ptr, edges, n_ternary, edge_ptr_ternary, edges_ternary);
  
  utils::binary_tree(root, n_ternary, edge_ptr_ternary, edges_ternary,
                     n_binary, child_ptr_binary, child_edges_binary, parent_binary);
}

// Build the RCTree for the rooted binary tree
template <class W>
void RCTree<W>::build(parlay::random_generator& gen, 
                      const parlay::sequence<W>& vertex_value,
                      const parlay::sequence<W>& edge_value) {
  edge_clusters.resize(n_binary - 1);
  parlay::parallel_for(0, n_binary, [&](uintV i) {
    if (i != root) {
      auto [p, id] = parent_binary[i];
      edge_clusters[id] = Cluster<W>(id, nullptr, nullptr, nullptr,
                                     nullptr, nullptr, p, i, id, id);
      if (id < n - 1)
        edge_clusters[id].set_val(edge_value[id]);
      else 
        edge_clusters[id].set_val(edge_base());
    }
  });
  
  vertex_clusters.resize(n_binary);
  parlay::parallel_for(0, n_binary, [&](uintV i) {
    vertex_clusters[i] = Cluster<W>(i, nullptr, nullptr, nullptr);
    if (i < n)
      vertex_clusters[i].set_val(vertex_value[i]);
    else
      vertex_clusters[i].set_val(vertex_base());
  });

  auto aux_ptr = parlay::sequence<Cluster<W>*>::from_function(
    n_binary - 1, [&](uintV i) {return &edge_clusters[i];});
  
  rc_clusters.resize(n_binary);

  auto deg = parlay::sequence<uintV>::from_function(
    n_binary, [&](uintV i) {return child_ptr_binary[i + 1] - child_ptr_binary[i] + (i != root);});

  parlay::sequence<uintV> state(n_binary), parent_peek(n_binary);

  auto rake = [&](uintV v) {
    Cluster<W> *binary_top = nullptr;
    Cluster<W> *unary[2] = {nullptr, nullptr};
    uintV usz = 0, boundary = 0, boundary_edge = 0;
    std::bernoulli_distribution dist(0.5);
    auto gen_v = gen[v];
    bool color_v = dist(gen_v);
    if (!color_v)  return;
    auto gen_parent = gen[parent_peek[v]];
    bool color_parent = dist(gen_parent) && state[parent_peek[v]];
    if (color_parent)  return;
    {
      auto [to, id] = parent_binary[v];
      auto &c = aux_ptr[id];
      binary_top = c;
      boundary = (c -> get_top_boundary()) ^ v ^ (c -> get_bottom_boundary());
      boundary_edge = (c -> get_top_boundary_edge()) ^ id ^ (c -> get_bottom_boundary_edge());
    }
    for (uintV i = child_ptr_binary[v]; i < child_ptr_binary[v + 1]; ++i)
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
    }
    aux_ptr[boundary_edge] = &rc_clusters[v];
    rc_clusters[v] = Cluster<W>(v, unary[0], unary[1], &vertex_clusters[v],
                                binary_top, boundary, boundary_edge);
    vertex_clusters[v].set_parent_cluster(&rc_clusters[v]);
    gbbs::write_add(deg.data() + rc_clusters[v].get_top_boundary(), -1);
    // deg[rc_clusters[v].get_top_boundary()]--;
    deg[v]--;
    rc_clusters[v].set_val(f_unary(&rc_clusters[v]));
  };

  auto compress = [&](uintV v) {
    Cluster<W> *binary[2] = {nullptr, nullptr}, *unary[2] = {nullptr, nullptr};
    uintV usz = 0, boundary[2] = {0, 0}, boundary_edge[2] = {0, 0};
    std::bernoulli_distribution dist(0.5);
    auto gen_v = gen[v];
    bool color_v = dist(gen_v);
    if (!color_v)  return;
    auto gen_parent = gen[parent_peek[v]];
    bool color_parent = dist(gen_parent) && state[parent_peek[v]];
    if (color_parent)  return;
    for (uintV i = child_ptr_binary[v]; i < child_ptr_binary[v + 1]; ++i)
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
    for (uintV i = child_ptr_binary[v]; i < child_ptr_binary[v + 1]; ++i)
    {
      auto [to, id] = child_edges_binary[i];
      auto &c = aux_ptr[id];
      c -> set_parent_cluster(&rc_clusters[v]);
      if (!(c -> is_binary()))
        unary[usz++] = c;
    }
    {
      auto [to, id] = parent_binary[v];
      auto &c = aux_ptr[id];
      c -> set_parent_cluster(&rc_clusters[v]);
    }
    aux_ptr[boundary_edge[0]] = &rc_clusters[v];
    aux_ptr[boundary_edge[1]] = &rc_clusters[v];
    rc_clusters[v] = Cluster<W>(v, unary[0], unary[1], &vertex_clusters[v],
                                binary[0], binary[1], boundary[0], boundary[1],
                                boundary_edge[0], boundary_edge[1]);
    vertex_clusters[v].set_parent_cluster(&rc_clusters[v]);
    deg[v] -= 2;
    rc_clusters[v].set_val(f_binary(&rc_clusters[v]));
  };
  
  while (true) {
    parlay::parallel_for(0, n_binary, [&](uintV i) {
      state[i] = 0;
      auto d = deg[i];
      if (i != root && d > 0 && d <= 2)
        state[i] = d;
    });
    
    if (parlay::all_of(state, [&](uintV i) {return i == 0;}))
      break;

    gen();
    parlay::parallel_for(0, n_binary, [&](uintV v) {
      if (state[v]) {
        auto [to, id] = parent_binary[v];
        auto &c = aux_ptr[id];
        parent_peek[v] = (c -> get_top_boundary()) ^ v ^ (c -> get_bottom_boundary());
      }
    });
    parlay::parallel_for(0, n_binary, [&](uintV v) {
      if (state[v] == 2)
        compress(v);
    });
    parlay::parallel_for(0, n_binary, [&](uintV v) {
      if (state[v] == 1)
        rake(v);
    });
  }

  auto finalize = [&](const uintV v) {
    Cluster<W> *unary[2] = {nullptr, nullptr};
    uintV usz = 0;
    for (uintV i = child_ptr_binary[v]; i < child_ptr_binary[v + 1]; ++i) {
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
  parlay::parallel_for(0, rc_clusters.size(), [&](const uintV i) {
    rc_clusters[i].walk();
  });
  parlay::parallel_for(0, vertex_clusters.size(), [&](const uintV i) {
    vertex_clusters[i].walk_from_parent_cluster();
  });
  parlay::parallel_for(0, edge_clusters.size(), [&](const uintV i) {
    edge_clusters[i].walk_from_parent_cluster();
  });
  // rc_clusters[root].print();
}

// Walking up the tree from leaf cluster to reevaluate value
template <class W>
void RCTree<W>::reevaluate(Cluster<W> *c) {
  while (c = c -> get_parent_cluster()) {
    switch(c -> get_cluster_type()) {
      case 0: c -> set_val(f_nullary(c)); break;
      case 1: c -> set_val(f_unary(c)); break;
      case 2: c -> set_val(f_binary(c)); break;
      default: assert(false);
    }
  }
}