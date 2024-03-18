#pragma once

#include <iostream>

#include "gbbs/macros.h"
#include "parlay/sequence.h"
#include "../utils/utils.h"

template <class T, class W>
class Cluster {
public:
  Cluster();
  Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster<T, W> *rep);
  Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster<T, W> *rep,
          Cluster<T, W> *bt, T b0, T be0);
  Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster *rep,
          Cluster<T, W> *bt, Cluster<T, W> *bb, T b0, T b1, T be0, T be1);
  bool is_binary() const;
  void print() const;
  void setval(const W& rval) {val = rval;}

  Cluster<T, W> *parent, *unary_cluster[2], *binary_cluster[2], *representative_cluster;
  T boundary[2], boundary_edge[2];
  // representative is the vertex/edge id if the cluster is an original vertex/edge.
  T representative;
  W val;
  int cluster_type;
};

// Default
template <class T, class W>
Cluster<T, W>::Cluster() {
  parent = nullptr;
  unary_cluster[0] = unary_cluster[1] = nullptr;
  binary_cluster[0] = binary_cluster[1] = nullptr;
  representative_cluster = nullptr;
  boundary[0] = boundary[1] = boundary_edge[0] = boundary_edge[1] = -1;
  representative = -1;
  cluster_type = -1;
}

// Nullary
template <class T, class W>
Cluster<T, W>::Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster<T, W> *rep) {
  parent = nullptr;
  unary_cluster[0] = u0, unary_cluster[1] = u1;
  binary_cluster[0] = binary_cluster[1] = nullptr;
  representative_cluster = rep;
  boundary[0] = boundary[1] = boundary_edge[0] = boundary_edge[1] = -1;
  representative = u;
  cluster_type = 0;
}

// Unary
template <class T, class W>
Cluster<T, W>::Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster<T, W> *rep,
                    Cluster<T, W> *bt, T b0, T be0) {
  parent = nullptr;
  unary_cluster[0] = u0, unary_cluster[1] = u1;
  binary_cluster[0] = bt, binary_cluster[1] = nullptr;
  representative_cluster = rep;
  boundary[0] = b0, boundary[1] = -1;
  boundary_edge[0] = be0, boundary_edge[1] = -1;
  representative = u;
  cluster_type = 1;
}

// Binary
template <class T, class W>
Cluster<T, W>::Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster *rep,
                    Cluster<T, W> *bt, Cluster<T, W> *bb, T b0, T b1, T be0, T be1) {
  parent = nullptr;
  unary_cluster[0] = u0, unary_cluster[1] = u1;
  binary_cluster[0] = bt, binary_cluster[1] = bb;
  representative_cluster = rep;
  boundary[0] = b0, boundary[1] = b1;
  boundary_edge[0] = be0, boundary_edge[1] = be1;
  representative = u;
  cluster_type = 2;
}

template <class T, class W>
bool Cluster<T, W>::is_binary() const {return cluster_type == 2;}

template <class T, class W>
void Cluster<T, W>::print() const {
  if (cluster_type == 0) {
    if (unary_cluster[0]) {
      std::cout << "root cluster vertex = " << representative << "\n";
      assert(this == representative_cluster -> parent);
      representative_cluster -> print();
      assert(this == unary_cluster[0] -> parent), unary_cluster[0] -> print();
      if (unary_cluster[1])
        assert(this == unary_cluster[1] -> parent), unary_cluster[1] -> print();
    } else {
      std::cout << "nullary cluster vertex_id = " << representative << "\n";
    }
  } else if (cluster_type == 1) {
    std::cout << "unary cluster vertex = " << representative << "\n";
    assert(this == representative_cluster -> parent);
    representative_cluster -> print();
    assert(this == binary_cluster[0] -> parent);
    binary_cluster[0] -> print();
    if (unary_cluster[0])
      assert(this == unary_cluster[0] -> parent), unary_cluster[0] -> print();
    if (unary_cluster[1])
      assert(this == unary_cluster[1] -> parent), unary_cluster[1] -> print();
  } else if (cluster_type == 2) {
    if (binary_cluster[0]) {
      std::cout << "binary cluster vertex = " << representative << "\n";
      assert(this == representative_cluster -> parent);
      representative_cluster -> print();
      assert(this == binary_cluster[0] -> parent);
      binary_cluster[0] -> print();
      assert(this == binary_cluster[1] -> parent);
      binary_cluster[1] -> print();
      if (unary_cluster[0])
        assert(this == unary_cluster[0] -> parent), unary_cluster[0] -> print();
      assert(unary_cluster[1] == nullptr);
    } else {
      std::cout << "binary cluster edge_id = " << representative << "\n";
    }
  } else {
    assert(false);
  }
}