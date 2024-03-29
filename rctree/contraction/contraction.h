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
  Contraction_Type<W> operator+(const Contraction_Type<W>& rhs) const {
    return Contraction_Type<W>(t - rhs.t, b - rhs.b, 0, 0, joined && rhs.joined);
  }
};

template <class W>
struct Contraction_MixOp {
  size_t timestamp;
  int type;
  uint u;
  W w;
  Contraction_MixOp(size_t timestamp = 0, int type = 0, uint u = 0, W w = 0) :
                    timestamp(timestamp), type(type), u(u), w(w) {}
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
  // void batch_operations(const parlay::sequence<Contraction_MixOp>& ops);
};

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_nullary(Cluster<T, Contraction_Type<W>> *c) {
  W s = c -> get_representative_cluster() -> get_val().t;
  for (int i = 0;i < 2; ++i)
    if (c -> get_unary_cluster(i))
      s += c -> get_unary_cluster(i) -> get_val().t;
  return Contraction_Type<W>(s, 0, 0, 0, true);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_unary(Cluster<T, Contraction_Type<W>> *c) {
  if (!c -> get_binary_top_cluster() -> get_val().joined)  
    return Contraction_Type<W>(c -> get_binary_top_cluster() -> get_val().t, 0, 
                               std::max(c -> get_binary_top_cluster() -> get_val().mt, c -> get_binary_top_cluster() -> get_val().mb), 0, true);
  W s = c -> get_representative_cluster() -> get_val().t;
  for (int i = 0;i < 2; ++i)
    if (c -> get_unary_cluster(i))
      s += c -> get_unary_cluster(i) -> get_val().t;
  return Contraction_Type<W>(s + c -> get_binary_top_cluster() -> get_val().t, 0, 
                             std::max(c -> get_binary_top_cluster() -> get_val().mt, c -> get_binary_top_cluster() -> get_val().mb), 0, true);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_binary(Cluster<T, Contraction_Type<W>> *c) {
  W s = c -> get_representative_cluster() -> get_val().t;
  for (int i = 0;i < 2; ++i)
    if (c -> get_unary_cluster(i))
      s += c -> get_unary_cluster(i) -> get_val().t;
  if (c -> get_binary_top_cluster() -> get_val().joined) {
    if (c -> get_binary_bottom_cluster() -> get_val().joined)
      return Contraction_Type<W>(s + c -> get_binary_top_cluster() -> get_val().t + c -> get_binary_bottom_cluster() -> get_val().t, 0, 
                                 std::max(c -> get_binary_top_cluster() -> get_val().mt, c -> get_binary_top_cluster() -> get_val().mb),
                                 std::max(c -> get_binary_bottom_cluster() -> get_val().mt, c -> get_binary_bottom_cluster() -> get_val().mb), true);
    return Contraction_Type<W>(c -> get_binary_top_cluster() -> get_val().t + c -> get_binary_bottom_cluster() -> get_val().t + s, 
                               c -> get_binary_bottom_cluster() -> get_val().b, 
                               std::max(c -> get_binary_top_cluster() -> get_val().mt, c -> get_binary_top_cluster() -> get_val().mb),
                               std::max(c -> get_binary_bottom_cluster() -> get_val().mt, c -> get_binary_bottom_cluster() -> get_val().mb), false);
  }
  if (c -> get_binary_bottom_cluster() -> get_val().joined)
    return Contraction_Type<W>(c -> get_binary_top_cluster() -> get_val().t, 
                               c -> get_binary_top_cluster() -> get_val().b + c -> get_binary_bottom_cluster() -> get_val().t + s, 
                               std::max(c -> get_binary_top_cluster() -> get_val().mt, c -> get_binary_top_cluster() -> get_val().mb),
                               std::max(c -> get_binary_bottom_cluster() -> get_val().mt, c -> get_binary_bottom_cluster() -> get_val().mb), false);
  return Contraction_Type<W>(c -> get_binary_top_cluster() -> get_val().t, c -> get_binary_bottom_cluster() -> get_val().b, 
                             std::max(c -> get_binary_top_cluster() -> get_val().mt, c -> get_binary_top_cluster() -> get_val().mb),
                             std::max(c -> get_binary_bottom_cluster() -> get_val().mt, c -> get_binary_bottom_cluster() -> get_val().mb), false);
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
  auto val = RCTree<T, Contraction_Type<W>>::vertex_clusters[u].get_val();
  val.t -= w;
  RCTree<T, Contraction_Type<W>>::vertex_clusters[u].set_val(val);
  RCTree<T, Contraction_Type<W>>::reevaluate(&RCTree<T, Contraction_Type<W>>::vertex_clusters[u]);
}

template <class T, class W>
void Contraction_RCTree<T, W>::join_edge(T e) {
  auto val = RCTree<T, Contraction_Type<W>>::edge_clusters[e].get_val();
  val.joined = true;
  RCTree<T, Contraction_Type<W>>::edge_clusters[e].set_val(val);
  RCTree<T, Contraction_Type<W>>::reevaluate(&RCTree<T, Contraction_Type<W>>::edge_clusters[e]);
}

template <class T, class W>
W Contraction_RCTree<T, W>::query_vertex_sum(T u) {
  auto c = RCTree<T, Contraction_Type<W>>::vertex_clusters[u].get_parent_cluster();
  while (true) {
    if (c -> get_binary_top_cluster() && c -> get_binary_top_cluster() -> get_val().joined) {
      c = RCTree<T, Contraction_Type<W>>::vertex_clusters[c -> get_top_boundary()].get_parent_cluster();
      continue;
    }
    if (c -> get_binary_bottom_cluster() && c -> get_binary_bottom_cluster() -> get_val().joined) {
      c = RCTree<T, Contraction_Type<W>>::vertex_clusters[c -> get_bottom_boundary()].get_parent_cluster();
      continue;
    }
    break;
  }
  W s = c -> get_representative_cluster() -> get_val().t;
  for (int i = 0; i < 2; ++i)
    if (c -> get_unary_cluster(i))
      s += c -> get_unary_cluster(i) -> get_val().t;
  if (c -> get_binary_top_cluster())
    s += c -> get_binary_top_cluster() -> get_val().b;
  if (c -> get_binary_bottom_cluster())
    s += c -> get_binary_bottom_cluster() -> get_val().t;
  return s;
}

template <class T, class W>
void Contraction_RCTree<T, W>::extract_max_path(Cluster<T, Contraction_Type<W>> *c, std::pair<uint, uint>& mx) {
  if (c -> get_parent_cluster() -> get_cluster_type() == 1) {
    if (c -> get_cluster_type() == 2)
      mx = {mx.first, 0};
    else
      mx = {std::max(mx.first, c -> get_parent_cluster() -> get_val().mt), 0};
  } else {
    if (c -> get_cluster_type() == 2) {
      if (c -> get_parent_cluster() -> get_binary_top_cluster() == c)
        mx = {mx.first, std::max(mx.second, c -> get_parent_cluster() -> get_val().mb)};
      else
        mx = {std::max(mx.first, c -> get_parent_cluster() -> get_val().mt), mx.second};
    } else {
      mx = {std::max(mx.first, c -> get_parent_cluster() -> get_val().mt), std::max(mx.first, c -> get_parent_cluster() -> get_val().mb)};
    }
  }
}

template <class T, class W>
uint Contraction_RCTree<T, W>::query_max_path(T u, T v) {
  uint du = 0, dv = 0;
  auto cu = &RCTree<T, Contraction_Type<W>>::vertex_clusters[u];
  auto cv = &RCTree<T, Contraction_Type<W>>::vertex_clusters[v];
  while (cu -> get_parent_cluster())  ++du, cu = cu -> get_parent_cluster();
  while (cv -> get_parent_cluster())  ++dv, cv = cv -> get_parent_cluster();
  cu = &RCTree<T, Contraction_Type<W>>::vertex_clusters[u];
  cv = &RCTree<T, Contraction_Type<W>>::vertex_clusters[v];
  if (du < dv)
    std::swap(cu, cv), std::swap(du, dv);
  std::pair<uint, uint> mu = {0, 0}, mv = {0, 0}; 
  while(du > dv)  --du, extract_max_path(cu, mu), cu = cu -> get_parent_cluster();

  while(cu -> get_parent_cluster() != cv -> get_parent_cluster()) {
    extract_max_path(cu, mu), extract_max_path(cv, mv);
    cu = cu -> get_parent_cluster(), cv = cv -> get_parent_cluster();
  }
  std::pair<uint, uint> bu, bv;
  if (cu -> get_cluster_type() == 0) {
    bu = {cu -> get_representative(), cu -> get_representative()};
  } else if (cu -> get_cluster_type() == 1) {
    bu = {cu -> get_top_boundary(), -1};
  } else {
    bu = {cu -> get_top_boundary(), cu -> get_bottom_boundary()};
  }
  if (cv -> get_cluster_type() == 0) {
    bv = {cv -> get_representative(), cv -> get_representative()};
  } else if (cv -> get_cluster_type() == 1) {
    bv = {cv -> get_top_boundary(), -1};
  } else {
    bv = {cv -> get_top_boundary(), cv -> get_bottom_boundary()};
  }
  if (bu.first == bv.first)
    return std::max(mu.first, mv.first);
  else if (bu.first == bv.second)
    return std::max(mu.first, mv.second);
  else if (bu.second == bv.first)
    return std::max(mu.second, mv.first);
  return std::max(mu.second, mv.second);
}

// template <class T, class W>
// void Contraction_RCTree<T, W>::batch_operations(const parlay::sequence<Contraction_MixOp>& ops) {

// }