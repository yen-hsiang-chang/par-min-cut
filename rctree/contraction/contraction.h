#pragma once
#include "../cluster.h"
#include "../rctree.h"

template <class T, class W>
class Contraction_RCTree : public RCTree<T, W> {
public:
  Contraction_RCTree(T n, const parlay::sequence<std::pair<T, T>>& edge_list)
    : RCTree<T, W>(n, edge_list) {}
  W f_nullary(const Cluster<T, W>& c) {return W();}
  W f_unary(const Cluster<T, W>& c) {return W();}
  W f_binary(const Cluster<T, W>& c) {return W();}
  W vertex_base() {return W();}
  W edge_base() {return W();}
};