#pragma once

#include <iostream>

#include "gbbs/macros.h"
#include "parlay/sequence.h"
#include "../utils/utils.h"

template <class T, class W>
class Cluster {
public:
  // Default constructor
  Cluster();
  // Nullary constructor
  Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster<T, W> *rep);
  // Unary constructor
  Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster<T, W> *rep,
          Cluster<T, W> *bt, T b0, T be0);
  // Binary constructor
  Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster *rep,
          Cluster<T, W> *bt, Cluster<T, W> *bb, T b0, T b1, T be0, T be1);
  // Check if this is a binary cluster
  bool is_binary() const;
  // Check the parent child relationship
  // 0: representative cluster
  // 1: binary top cluster, 2: binary bottom cluster
  // 3: first unary cluster, 4: second unary cluster
  int parent_child_relation() const;
  // Get top/bottom boundary
  T get_top_boundary() const;
  T get_bottom_boundary() const;
  // Get top/bottom boundary edge
  T get_top_boundary_edge() const;
  T get_bottom_boundary_edge() const;
  // Get representative
  T get_representative() const;
  // Get cluster type
  int get_cluster_type() const;
  // Set parent cluster
  void set_parent_cluster(Cluster<T, W> *p);
  // Get parent cluster
  Cluster<T, W>* get_parent_cluster() const;
  // Get representative cluster
  Cluster<T, W>* get_representative_cluster() const;
  // Get binary top cluster
  Cluster<T, W>* get_binary_top_cluster() const;
  // Get binary bottom cluster
  Cluster<T, W>* get_binary_bottom_cluster() const;
  // Get unary cluster by idx
  Cluster<T, W>* get_unary_cluster(T idx) const;
  // Get maintained value
  W get_val() const;
  // Set maintained value
  void set_val(const W& rval);
  // Get level
  int get_level() const;
  // Get tid
  T get_tid() const;
  // Walking up the tree to set level and tid
  void walk();
  // Set level and tid from parent, only for leaf clusters
  void walk_from_parent_cluster();
  // Print debug info for the rctree structure through dfs
  void print() const;

private:
  // Pointer to parent cluster, binary clusters, unary clusters and representative cluster
  Cluster<T, W> *parent_cluster;
  Cluster<T, W> *binary_cluster[2];
  Cluster<T, W> *unary_cluster[2];
  Cluster<T, W> *representative_cluster;
  // Boundary vertices and edges
  T boundary[2], boundary_edge[2];
  // Representative is the vertex/edge id if the cluster is an original vertex/edge
  // Otherwise, it is the id of the representative vertex
  T representative;
  // Maintained value
  W val;
  // Cluster type. 0: Nullary, 1: Unary, 2: Binary
  int cluster_type;
  // Level in the rctree
  int level;
  // Number of clusters in the subtree
  T size;
  // Timestamp in preorder traversal
  T tid;
};

// Default constructor
template <class T, class W>
Cluster<T, W>::Cluster() {
  parent_cluster = nullptr;
  unary_cluster[0] = unary_cluster[1] = nullptr;
  binary_cluster[0] = binary_cluster[1] = nullptr;
  representative_cluster = nullptr;
  boundary[0] = boundary[1] = boundary_edge[0] = boundary_edge[1] = -1;
  representative = -1;
  cluster_type = -1;
  size = 0;
}

// Nullary constructor
template <class T, class W>
Cluster<T, W>::Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster<T, W> *rep) {
  parent_cluster = nullptr;
  unary_cluster[0] = u0, unary_cluster[1] = u1;
  binary_cluster[0] = binary_cluster[1] = nullptr;
  representative_cluster = rep;
  boundary[0] = boundary[1] = boundary_edge[0] = boundary_edge[1] = -1;
  representative = u;
  cluster_type = 0;
  size = 1 + (rep ? rep -> size : 0) + (u0 ? u0 -> size : 0) + (u1 ? u1 -> size : 0);
}

// Unary constructor
template <class T, class W>
Cluster<T, W>::Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster<T, W> *rep,
                    Cluster<T, W> *bt, T b0, T be0) {
  parent_cluster = nullptr;
  unary_cluster[0] = u0, unary_cluster[1] = u1;
  binary_cluster[0] = bt, binary_cluster[1] = nullptr;
  representative_cluster = rep;
  boundary[0] = b0, boundary[1] = -1;
  boundary_edge[0] = be0, boundary_edge[1] = -1;
  representative = u;
  cluster_type = 1;
  size = 1 + (rep ? rep -> size : 0) + (u0 ? u0 -> size : 0) + (u1 ? u1 -> size : 0) + (bt ? bt -> size : 0);
}

// Binary constructor
template <class T, class W>
Cluster<T, W>::Cluster(T u, Cluster<T, W> *u0, Cluster<T, W> *u1, Cluster *rep,
                    Cluster<T, W> *bt, Cluster<T, W> *bb, T b0, T b1, T be0, T be1) {
  parent_cluster = nullptr;
  unary_cluster[0] = u0, unary_cluster[1] = u1;
  binary_cluster[0] = bt, binary_cluster[1] = bb;
  representative_cluster = rep;
  boundary[0] = b0, boundary[1] = b1;
  boundary_edge[0] = be0, boundary_edge[1] = be1;
  representative = u;
  cluster_type = 2;
  size = 1 + (rep ? rep -> size : 0) + (u0 ? u0 -> size : 0) + (u1 ? u1 -> size : 0) + (bt ? bt -> size : 0) + (bb ? bb -> size : 0);
}

// Check if this is a binary cluster
template <class T, class W>
inline bool Cluster<T, W>::is_binary() const {
  return cluster_type == 2;
}

// Check the parent child relationship
// 0: representative cluster
// 1: binary top cluster, 2: binary bottom cluster
// 3: first unary cluster, 4: second unary cluster
template <class T, class W>
int Cluster<T, W>::parent_child_relation() const {
  if (cluster_type == 0)  return 0;
  if (cluster_type == 2) {
    if (this == parent_cluster -> binary_cluster[0])
      return 1;
    return 2;
  }
  if (this == parent_cluster -> unary_cluster[0])
    return 3;
  return 4;
}

// Get top boundary
template <class T, class W>
inline T Cluster<T, W>::get_top_boundary() const {
  return boundary[0];
}

// Get bottom boundary
template <class T, class W>
inline T Cluster<T, W>::get_bottom_boundary() const {
  return boundary[1];
}

// Get top boundary edge
template <class T, class W>
inline T Cluster<T, W>::get_top_boundary_edge() const {
  return boundary_edge[0];
}

// Get bottom boundary edge
template <class T, class W>
inline T Cluster<T, W>::get_bottom_boundary_edge() const {
  return boundary_edge[1];
}

// Get representative
template <class T, class W>
inline T Cluster<T, W>::get_representative() const {
  return representative;
}

// Get cluster type
template <class T, class W>
inline int Cluster<T, W>::get_cluster_type() const {
  return cluster_type;
}

// Set parent cluster
template <class T, class W>
inline void Cluster<T, W>::set_parent_cluster(Cluster<T, W> *p) {
  parent_cluster = p;
}

// Get parent cluster
template <class T, class W>
inline Cluster<T, W>* Cluster<T, W>::get_parent_cluster() const {
  return parent_cluster;
}

// Get representative cluster
template <class T, class W>
inline Cluster<T, W>* Cluster<T, W>::get_representative_cluster() const {
  return representative_cluster;
}

// Get binary top cluster
template <class T, class W>
inline Cluster<T, W>* Cluster<T, W>::get_binary_top_cluster() const {
  return binary_cluster[0];
}

// Get binary bottom cluster
template <class T, class W>
inline Cluster<T, W>* Cluster<T, W>::get_binary_bottom_cluster() const {
  return binary_cluster[1];
}

// Get unary cluster by idx
template <class T, class W>
inline Cluster<T, W>* Cluster<T, W>::get_unary_cluster(T idx) const {
  return unary_cluster[idx];
}

// Get maintained value
template <class T, class W>
inline W Cluster<T, W>::get_val() const {
  return val;
}

// Set maintained value
template <class T, class W>
inline void Cluster<T, W>::set_val(const W& rval) {
  val = rval;
}

// Get level
template <class T, class W>
inline int Cluster<T, W>::get_level() const {
  return level;
}

// Get tid
template <class T, class W>
inline T Cluster<T, W>::get_tid() const {
  return tid;
}

// Walking up the tree to set level and tid
template <class T, class W>
void Cluster<T, W>::walk() {
  auto c = this;
  level = 0;
  tid = 0;
  while (c -> parent_cluster) {
    auto nc = c -> parent_cluster;
    ++level;
    ++tid;
    int pc = c -> parent_child_relation();
    if (pc > 0) tid += nc -> representative_cluster -> size;
    if (pc > 1 && nc -> binary_cluster[0])  tid += nc -> binary_cluster[0] -> size;
    if (pc > 2 && nc -> binary_cluster[1])  tid += nc -> binary_cluster[1] -> size;
    if (pc > 3 && nc -> unary_cluster[0]) tid += nc -> unary_cluster[0] -> size;
    c = nc;
  }
}

// Set level and tid from parent, only for leaf clusters
template <class T, class W>
inline void Cluster<T, W>::walk_from_parent_cluster() {
  level = parent_cluster -> level + 1;
  if (cluster_type == 0) {
    tid = parent_cluster -> tid + 1;
  } else {
    if (parent_child_relation() == 1)
      tid = parent_cluster -> tid + 2;
    else
      tid = parent_cluster -> tid + 3;
  }
}

// Print debug info for the rctree structure through dfs
template <class T, class W>
void Cluster<T, W>::print() const {
  std::cout << "level = " << level << ", tid = " << tid << "\n";
  if (cluster_type == 0) {
    if (unary_cluster[0]) {
      std::cout << "root cluster vertex = " << representative << "\n";
      assert(this == representative_cluster -> parent_cluster);
      representative_cluster -> print();
      assert(this == unary_cluster[0] -> parent_cluster), unary_cluster[0] -> print();
      if (unary_cluster[1])
        assert(this == unary_cluster[1] -> parent_cluster), unary_cluster[1] -> print();
    } else {
      std::cout << "nullary cluster vertex_id = " << representative << "\n";
    }
  } else if (cluster_type == 1) {
    std::cout << "unary cluster vertex = " << representative << "\n";
    assert(this == representative_cluster -> parent_cluster);
    representative_cluster -> print();
    assert(this == binary_cluster[0] -> parent_cluster);
    binary_cluster[0] -> print();
    if (unary_cluster[0])
      assert(this == unary_cluster[0] -> parent_cluster), unary_cluster[0] -> print();
    if (unary_cluster[1])
      assert(this == unary_cluster[1] -> parent_cluster), unary_cluster[1] -> print();
  } else if (cluster_type == 2) {
    if (binary_cluster[0]) {
      std::cout << "binary cluster vertex = " << representative << "\n";
      assert(this == representative_cluster -> parent_cluster);
      representative_cluster -> print();
      assert(this == binary_cluster[0] -> parent_cluster);
      binary_cluster[0] -> print();
      assert(this == binary_cluster[1] -> parent_cluster);
      binary_cluster[1] -> print();
      if (unary_cluster[0])
        assert(this == unary_cluster[0] -> parent_cluster), unary_cluster[0] -> print();
      assert(unary_cluster[1] == nullptr);
    } else {
      std::cout << "binary cluster edge_id = " << representative << "\n";
    }
  } else {
    assert(false);
  }
}