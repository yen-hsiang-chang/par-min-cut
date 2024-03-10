#pragma once

#include <iostream>

#include "gbbs/macros.h"
#include "parlay/sequence.h"
#include "../utils/utils.h"

class Cluster {
public:
  Cluster();
  Cluster(uint u, const parlay::sequence<Cluster*>& un, Cluster *rep);
  Cluster(uint u, const parlay::sequence<Cluster*>& un, Cluster *rep,
          Cluster *bt, uint b0, uint be0);
  Cluster(uint u, const parlay::sequence<Cluster*>& un, Cluster *rep,
          Cluster *bt, Cluster *bb, uint b0, uint b1, uint be0, uint be1);
  bool is_binary() const;
  void print() const;

  Cluster *parent, *binary_cluster[2], *representative_cluster;
  parlay::sequence<Cluster*> unary;
  uint boundary[2], boundary_edge[2];
  // representative is the vertex/edge id if the cluster is an original vertex/edge.
  uint representative;
  int cluster_type;
};

// Default
Cluster::Cluster() {
  parent = binary_cluster[0] = binary_cluster[1] = nullptr;
  representative_cluster = nullptr;
  boundary[0] = boundary[1] = boundary_edge[0] = boundary_edge[1] = -1;
  representative = -1;
  cluster_type = -1;
}

// Nullary
Cluster::Cluster(uint u, const parlay::sequence<Cluster*>& un, Cluster *rep) {
  parent = binary_cluster[0] = binary_cluster[1] = nullptr;
  representative_cluster = rep;
  unary = un;
  boundary[0] = boundary[1] = boundary_edge[0] = boundary_edge[1] = -1;
  representative = u;
  cluster_type = 0;
}

// Unary
Cluster::Cluster(uint u, const parlay::sequence<Cluster*>& un, Cluster *rep,
                 Cluster *bt, uint b0, uint be0) {
  parent = binary_cluster[1] = nullptr;
  binary_cluster[0] = bt;
  representative_cluster = rep;
  unary = un;
  boundary[0] = b0;
  boundary[1] = -1;
  boundary_edge[0] = be0;
  boundary_edge[1] = -1;
  representative = u;
  cluster_type = 1;
}

// Binary
Cluster::Cluster(uint u, const parlay::sequence<Cluster*>& un, Cluster *rep,
                 Cluster *bt, Cluster *bb, uint b0, uint b1, uint be0, uint be1) {
  parent = nullptr;
  binary_cluster[0] = bt;
  binary_cluster[1] = bb;
  representative_cluster = rep;
  unary = un;
  boundary[0] = b0;
  boundary[1] = b1;
  boundary_edge[0] = be0;
  boundary_edge[1] = be1;
  representative = u;
  cluster_type = 2;
}

bool Cluster::is_binary() const {return cluster_type == 2;}

void Cluster::print() const {
  if (cluster_type == 0) {
    if (unary.size()) {
      std::cout << "root cluster vertex = " << representative << "\n";
      assert(this == representative_cluster -> parent);
      representative_cluster -> print();
      for (auto p : unary)
        assert(this == p -> parent), p -> print();
    } else {
      assert(parent != nullptr);
      std::cout << "nullary cluster vertex_id = " << representative << "\n";
    }
  } else if (cluster_type == 1) {
    assert(parent != nullptr);
    std::cout << "unary cluster vertex = " << representative << "\n";
    assert(this == representative_cluster -> parent);
    representative_cluster -> print();
    assert(this == binary_cluster[0] -> parent);
    binary_cluster[0] -> print();
    for (auto p : unary)
      assert(this == p -> parent), p -> print();
  } else if (cluster_type == 2) {
    assert(parent != nullptr);
    if (binary_cluster[0] != nullptr) {
      std::cout << "binary cluster vertex = " << representative << "\n";
      assert(this == representative_cluster -> parent);
      representative_cluster -> print();
      assert(this == binary_cluster[0] -> parent);
      binary_cluster[0] -> print();
      assert(this == binary_cluster[1] -> parent);
      binary_cluster[1] -> print();
      for (auto p : unary)
        assert(this == p -> parent), p -> print();
    } else {
      std::cout << "binary cluster edge_id = " << representative << "\n";
    }
  } else {
    assert(false);
  }
}