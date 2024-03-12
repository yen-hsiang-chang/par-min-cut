#pragma once
#include "../cluster.h"
#include "../rctree.h"

template <class W>
struct Contraction_Type {
  W t, b;
  bool joined;
  Contraction_Type(W t = 0, W b = 0, bool joined = false) :
                   t(t), b(b), joined(joined) {}
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
  void subtract_vertex_weight(T i, W w);
  void join_edge(T i);
  W query_vertex(T i);
};

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_nullary(Cluster<T, Contraction_Type<W>> *c) {
  W s = c -> representative_cluster -> val.t;
  for (int i = 0;i < 2; ++i)
    if (c -> unary_cluster[i] != nullptr)
      s += c -> unary_cluster[i] -> val.t;
  return Contraction_Type<W>(s, 0, true);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_unary(Cluster<T, Contraction_Type<W>> *c) {
  if (!c -> binary_cluster[0] -> val.joined)  return Contraction_Type<W>(c -> binary_cluster[0] -> val.t, 0, true);
  W s = c -> representative_cluster -> val.t;
  for (int i = 0;i < 2; ++i)
    if (c -> unary_cluster[i] != nullptr)
      s += c -> unary_cluster[i] -> val.t;
  return Contraction_Type<W>(s + c -> binary_cluster[0] -> val.t, 0, true);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_binary(Cluster<T, Contraction_Type<W>> *c) {
  W s = c -> representative_cluster -> val.t;
  for (int i = 0;i < 2; ++i)
    if (c -> unary_cluster[i] != nullptr)
      s += c -> unary_cluster[i] -> val.t;
  if (c -> binary_cluster[0] -> val.joined) {
    if (c -> binary_cluster[1] -> val.joined)
      return Contraction_Type<W>(s + c -> binary_cluster[0] -> val.t + c -> binary_cluster[1] -> val.t, 0, true);
    return Contraction_Type<W>(c -> binary_cluster[0] -> val.t + c -> binary_cluster[1] -> val.t + s, c -> binary_cluster[1] -> val.b, false);
  }
  if (c -> binary_cluster[1] -> val.joined)
    return Contraction_Type<W>(c -> binary_cluster[0] -> val.t, c -> binary_cluster[0] -> val.b + c -> binary_cluster[1] -> val.b+ s);
  return Contraction_Type<W>(c -> binary_cluster[0] -> val.t, c -> binary_cluster[1] -> val.b, false);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::vertex_base() {
  return Contraction_Type<W>(0, 0, true);
}

template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::edge_base() {
  return Contraction_Type<W>(0, 0, true);
}

template <class T, class W>
void Contraction_RCTree<T, W>::subtract_vertex_weight(T i, W w) {
  RCTree<T, Contraction_Type<W>>::vertex_clusters[i].val.t -= w;
  reevaluate(&RCTree<T, Contraction_Type<W>>::vertex_clusters[i]);
}

template <class T, class W>
void Contraction_RCTree<T, W>::join_edge(T i) {
  RCTree<T, Contraction_Type<W>>::edge_clusters[i].joined = true;
}

template <class T, class W>
W Contraction_RCTree<T, W>::query_vertex(T i) {
  auto c = RCTree<T, Contraction_Type<W>>::vertex_clusters[i].parent;
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
    if (c -> unary_cluster[i] != nullptr)
      s += c -> unary_cluster[i] -> val.t;
  if (c -> binary_cluster[0])
    s += c -> binary_cluster[0] -> val.b;
  if (c -> binary_cluster[1])
    s += c -> binary_cluster[1] -> val.t;
  return s;
}