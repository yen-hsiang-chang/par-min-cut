#pragma once
#include "../cluster.h"
#include "../rctree.h"

template <class W>
struct Contraction_Type {
  W t, b; // sum of weights of connected componenent from top/bottom boundaries, only count for current cluster 
  uint mt, mb; // max edge id from representative vertex to top/bottom boundaries
  bool joined;
  Contraction_Type(W t = 0, W b = 0, uint mt = 0, uint mb = 0, bool joined = false) :
                   t(t), b(b), mt(mt), mb(mb), joined(joined) {}
};

template <class T, class W>
class Contraction_RCTree : public RCTree<T, Contraction_Type<W>> {
public:
  Contraction_RCTree(T n, const parlay::sequence<std::pair<T, T>>& edge_list)
    : RCTree<T, Contraction_Type<W>>(n, edge_list) {}
  Contraction_Type<W> f_nullary(Cluster<T, Contraction_Type<W>> *c);
  Contraction_Type<W> f_unary(Cluster<T, Contraction_Type<W>> *c);
  Contraction_Type<W> f_binary(Cluster<T, Contraction_Type<W>> *c);
  Contraction_Type<W> vertex_base();
  Contraction_Type<W> edge_base();
  void subtract_vertex_weight(T u, W w);
  void join_edge(T e);
  W query_vertex_sum(T u);
  void extract_max_path(Cluster<T, Contraction_Type<W>> *c, std::pair<uint, uint>& mx);
  uint query_max_path(T u, T v);
};

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_nullary(Cluster<T, Contraction_Type<W>> *c) {
  W s = c -> representative_cluster -> val.t;
  for (int i = 0;i < 2; ++i)
    if (c -> unary_cluster[i])
      s += c -> unary_cluster[i] -> val.t;
  return Contraction_Type<W>(s, 0, 0, 0, true);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_unary(Cluster<T, Contraction_Type<W>> *c) {
  if (!c -> binary_cluster[0] -> val.joined)  
    return Contraction_Type<W>(c -> binary_cluster[0] -> val.t, 0, 
                               std::max(c -> binary_cluster[0] -> val.mt, c -> binary_cluster[0] -> val.mb), 0, true);
  W s = c -> representative_cluster -> val.t;
  for (int i = 0;i < 2; ++i)
    if (c -> unary_cluster[i])
      s += c -> unary_cluster[i] -> val.t;
  return Contraction_Type<W>(s + c -> binary_cluster[0] -> val.t, 0, 
                             std::max(c -> binary_cluster[0] -> val.mt, c -> binary_cluster[0] -> val.mb), 0, true);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_binary(Cluster<T, Contraction_Type<W>> *c) {
  W s = c -> representative_cluster -> val.t;
  for (int i = 0;i < 2; ++i)
    if (c -> unary_cluster[i])
      s += c -> unary_cluster[i] -> val.t;
  if (c -> binary_cluster[0] -> val.joined) {
    if (c -> binary_cluster[1] -> val.joined)
      return Contraction_Type<W>(s + c -> binary_cluster[0] -> val.t + c -> binary_cluster[1] -> val.t, 0, 
                                 std::max(c -> binary_cluster[0] -> val.mt, c -> binary_cluster[0] -> val.mb),
                                 std::max(c -> binary_cluster[1] -> val.mt, c -> binary_cluster[1] -> val.mb), true);
    return Contraction_Type<W>(c -> binary_cluster[0] -> val.t + c -> binary_cluster[1] -> val.t + s, 
                               c -> binary_cluster[1] -> val.b, 
                               std::max(c -> binary_cluster[0] -> val.mt, c -> binary_cluster[0] -> val.mb),
                               std::max(c -> binary_cluster[1] -> val.mt, c -> binary_cluster[1] -> val.mb), false);
  }
  if (c -> binary_cluster[1] -> val.joined)
    return Contraction_Type<W>(c -> binary_cluster[0] -> val.t, 
                               c -> binary_cluster[0] -> val.b + c -> binary_cluster[1] -> val.b + s, 
                               std::max(c -> binary_cluster[0] -> val.mt, c -> binary_cluster[0] -> val.mb),
                               std::max(c -> binary_cluster[1] -> val.mt, c -> binary_cluster[1] -> val.mb), false);
  return Contraction_Type<W>(c -> binary_cluster[0] -> val.t, c -> binary_cluster[1] -> val.b, 
                             std::max(c -> binary_cluster[0] -> val.mt, c -> binary_cluster[0] -> val.mb),
                             std::max(c -> binary_cluster[1] -> val.mt, c -> binary_cluster[1] -> val.mb), false);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::vertex_base() {
  return Contraction_Type<W>(0, 0, 0, 0, true);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::edge_base() {
  return Contraction_Type<W>(0, 0, 0, 0, true);
}

template <class T, class W>
void Contraction_RCTree<T, W>::subtract_vertex_weight(T u, W w) {
  RCTree<T, Contraction_Type<W>>::vertex_clusters[u].val.t -= w;
  reevaluate(&RCTree<T, Contraction_Type<W>>::vertex_clusters[u]);
}

template <class T, class W>
void Contraction_RCTree<T, W>::join_edge(T e) {
  RCTree<T, Contraction_Type<W>>::edge_clusters[e].joined = true;
  reevaluate(&RCTree<T, Contraction_Type<W>>::edge_clusters[e]);
}

template <class T, class W>
W Contraction_RCTree<T, W>::query_vertex_sum(T u) {
  auto c = RCTree<T, Contraction_Type<W>>::vertex_clusters[u].parent;
  while (true) {
    if (c -> binary_cluster[0] && c -> binary_cluster[0] -> val.joined) {
      c = RCTree<T, Contraction_Type<W>>::vertex_clusters[c -> boundary[0]].parent;
      continue;
    }
    if (c -> binary_cluster[1] && c -> binary_cluster[1] -> val.joined) {
      c = RCTree<T, Contraction_Type<W>>::vertex_clusters[c -> boundary[1]].parent;
      continue;
    }
    break;
  }
  W s = c -> representative_cluster -> val.t;
  for (int i = 0;i < 2; ++i)
    if (c -> unary_cluster[i])
      s += c -> unary_cluster[i] -> val.t;
  if (c -> binary_cluster[0])
    s += c -> binary_cluster[0] -> val.b;
  if (c -> binary_cluster[1])
    s += c -> binary_cluster[1] -> val.t;
  return s;
}

template <class T, class W>
void Contraction_RCTree<T, W>::extract_max_path(Cluster<T, Contraction_Type<W>> *c, std::pair<uint, uint>& mx) {
  if (c -> parent -> cluster_type == 1) {
    if (c -> cluster_type == 2)
      mx = {mx.first, 0};
    else
      mx = {std::max(mx.first, c -> parent -> val.mt), 0};
  } else {
    if (c -> cluster_type == 2) {
      if (c -> parent -> binary_cluster[0] == c)
        mx = {mx.first, std::max(mx.second, c -> parent -> val.mb)};
      else
        mx = {std::max(mx.first, c -> parent -> val.mt), mx.second};
    } else {
      mx = {std::max(mx.first, c -> parent -> val.mt), std::max(mx.first, c -> parent -> val.mb)};
    }
  }
}

template <class T, class W>
uint Contraction_RCTree<T, W>::query_max_path(T u, T v) {
  uint du = 0, dv = 0;
  auto cu = &RCTree<T, Contraction_Type<W>>::vertex_clusters[u];
  auto cv = &RCTree<T, Contraction_Type<W>>::vertex_clusters[v];
  while (cu -> parent)  ++du, cu = cu -> parent;
  while (cv -> parent)  ++dv, cv = cv -> parent;
  cu = &RCTree<T, Contraction_Type<W>>::vertex_clusters[u];
  cv = &RCTree<T, Contraction_Type<W>>::vertex_clusters[v];
  if (du < dv)
    std::swap(cu, cv), std::swap(du, dv);
  std::pair<uint, uint> mu = {0, 0}, mv = {0, 0}; 
  while(du > dv)  --du, extract_max_path(cu, mu), cu = cu -> parent;

  while(cu -> parent != cv -> parent) {
    extract_max_path(cu, mu), extract_max_path(cv, mv);
    cu = cu -> parent, cv = cv -> parent;
  }
  std::pair<uint, uint> bu, bv;
  if (cu -> cluster_type == 0) {
    bu = {cu -> representative, cu -> representative};
  } else if (cu -> cluster_type == 1) {
    bu = {cu -> boundary[0], -1};
  } else {
    bu = {cu -> boundary[0], cu -> boundary[1]};
  }
  if (cv -> cluster_type == 0) {
    bv = {cv -> representative, cv -> representative};
  } else if (cv -> cluster_type == 1) {
    bv = {cv -> boundary[0], -1};
  } else {
    bv = {cv -> boundary[0], cv -> boundary[1]};
  }
  if (bu.first == bv.first)
    return std::max(mu.first, mv.first);
  else if (bu.first == bv.second)
    return std::max(mu.first, mv.second);
  else if (bu.second == bv.first)
    return std::max(mu.second, mv.first);
  return std::max(mu.second, mv.second);
}