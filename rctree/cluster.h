#pragma once

#include <iostream>

#include "gbbs/macros.h"
#include "parlay/sequence.h"
#include "../utils/utils.h"

template <class W>
class Cluster {
public:
  // Default constructor
  Cluster();
  // Nullary constructor
  Cluster(uintV u, Cluster<W> *u0, Cluster<W> *u1, Cluster<W> *rep);
  // Unary constructor
  Cluster(uintV u, Cluster<W> *u0, Cluster<W> *u1, Cluster<W> *rep,
          Cluster<W> *bt, uintV b0, uintV be0);
  // Binary constructor
  Cluster(uintV u, Cluster<W> *u0, Cluster<W> *u1, Cluster *rep,
          Cluster<W> *bt, Cluster<W> *bb, uintV b0, uintV b1, uintV be0, uintV be1);
  // Check if this is a binary cluster
  bool is_binary() const;
  // Check the parent child relationship by tid
  // 0: representative cluster
  // 1: binary top cluster, 2: binary bottom cluster
  // 3: first unary cluster, 4: second unary cluster
  int parent_child_relation_by_tid() const;
  // Get top/bottom boundary
  uintV get_top_boundary() const;
  uintV get_bottom_boundary() const;
  // Get top/bottom boundary edge
  uintV get_top_boundary_edge() const;
  uintV get_bottom_boundary_edge() const;
  // Get representative
  uintV get_representative() const;
  // Get cluster type
  int get_cluster_type() const;
  // Set parent cluster
  void set_parent_cluster(Cluster<W> *c);
  // Get parent cluster
  Cluster<W>* get_parent_cluster() const;
  // Set representative cluster
  void set_representative_cluster(Cluster<W> *c);
  // Get representative cluster
  Cluster<W>* get_representative_cluster() const;
  // Set binary top cluster
  void set_binary_top_cluster(Cluster<W> *c);
  // Get binary top cluster
  Cluster<W>* get_binary_top_cluster() const;
  // Set binary bottom cluster
  void set_binary_bottom_cluster(Cluster<W> *c);
  // Get binary bottom cluster
  Cluster<W>* get_binary_bottom_cluster() const;
  // Set unary cluster by idx
  void set_unary_cluster(int idx, Cluster<W> *c);
  // Get unary cluster by idx
  Cluster<W>* get_unary_cluster(int idx) const;
  // Get maintained value
  W get_val() const;
  // Set maintained value
  void set_val(const W& rval);
  // Get level
  int get_level() const;
  // Get tid
  uint64_t get_tid() const;
  // Walking up the tree to set level and tid
  void walk();
  // Set level and tid from parent, only for leaf clusters
  void walk_from_parent_cluster();
  // Print debug info for the rctree structure through dfs
  void print() const;

private:
  // Pointer to parent cluster, binary clusters, unary clusters and representative cluster
  Cluster<W> *parent_cluster;
  Cluster<W> *binary_cluster[2];
  Cluster<W> *unary_cluster[2];
  Cluster<W> *representative_cluster;
  // Boundary vertices and edges
  uintV boundary[2], boundary_edge[2];
  // Representative is the vertex/edge id if the cluster is an original vertex/edge
  // Otherwise, it is the id of the representative vertex
  uintV representative;
  // Maintained value
  W val;
  // Cluster type. 0: Nullary, 1: Unary, 2: Binary
  int cluster_type;
  // Level in the rctree
  int level;
  // Number of clusters in the subtree
  uint64_t size;
  // Timestamp in preorder traversal
  uint64_t tid;
  // Check the parent child relationship by pointer
  // 0: representative cluster
  // 1: binary top cluster, 2: binary bottom cluster
  // 3: first unary cluster, 4: second unary cluster
  int parent_child_relation_by_pointer() const;
};

// Default constructor
template <class W>
Cluster<W>::Cluster() {
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
template <class W>
Cluster<W>::Cluster(uintV u, Cluster<W> *u0, Cluster<W> *u1, Cluster<W> *rep) {
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
template <class W>
Cluster<W>::Cluster(uintV u, Cluster<W> *u0, Cluster<W> *u1, Cluster<W> *rep,
                    Cluster<W> *bt, uintV b0, uintV be0) {
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
template <class W>
Cluster<W>::Cluster(uintV u, Cluster<W> *u0, Cluster<W> *u1, Cluster *rep,
                    Cluster<W> *bt, Cluster<W> *bb, uintV b0, uintV b1, uintV be0, uintV be1) {
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
template <class W>
inline bool Cluster<W>::is_binary() const {
  return cluster_type == 2;
}

// Check the parent child relationship by tid
// 0: representative cluster
// 1: binary top cluster, 2: binary bottom cluster
// 3: first unary cluster, 4: second unary cluster
template <class W>
int Cluster<W>::parent_child_relation_by_tid() const {
  if (cluster_type == 0)  return 0;
  if (cluster_type == 2) {
    if (tid == parent_cluster -> binary_cluster[0] -> tid)
      return 1;
    return 2;
  }
  if (tid == parent_cluster -> unary_cluster[0] -> tid)
    return 3;
  return 4;
}

// Check the parent child relationship by pointer
// 0: representative cluster
// 1: binary top cluster, 2: binary bottom cluster
// 3: first unary cluster, 4: second unary cluster
template <class W>
int Cluster<W>::parent_child_relation_by_pointer() const {
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
template <class W>
inline uintV Cluster<W>::get_top_boundary() const {
  return boundary[0];
}

// Get bottom boundary
template <class W>
inline uintV Cluster<W>::get_bottom_boundary() const {
  return boundary[1];
}

// Get top boundary edge
template <class W>
inline uintV Cluster<W>::get_top_boundary_edge() const {
  return boundary_edge[0];
}

// Get bottom boundary edge
template <class W>
inline uintV Cluster<W>::get_bottom_boundary_edge() const {
  return boundary_edge[1];
}

// Get representative
template <class W>
inline uintV Cluster<W>::get_representative() const {
  return representative;
}

// Get cluster type
template <class W>
inline int Cluster<W>::get_cluster_type() const {
  return cluster_type;
}

// Set parent cluster
template <class W>
inline void Cluster<W>::set_parent_cluster(Cluster<W> *c) {
  parent_cluster = c;
}

// Get parent cluster
template <class W>
inline Cluster<W>* Cluster<W>::get_parent_cluster() const {
  return parent_cluster;
}

// Set representative cluster
template <class W>
inline void Cluster<W>::set_representative_cluster(Cluster<W> *c) {
  representative_cluster = c;
}

// Get representative cluster
template <class W>
inline Cluster<W>* Cluster<W>::get_representative_cluster() const {
  return representative_cluster;
}

// Set binary top cluster
template <class W>
inline void Cluster<W>::set_binary_top_cluster(Cluster<W> *c) {
  binary_cluster[0] = c;
}

// Get binary top cluster
template <class W>
inline Cluster<W>* Cluster<W>::get_binary_top_cluster() const {
  return binary_cluster[0];
}

// Set binary bottom cluster
template <class W>
inline void Cluster<W>::set_binary_bottom_cluster(Cluster<W> *c) {
  binary_cluster[1] = c;
}

// Get binary bottom cluster
template <class W>
inline Cluster<W>* Cluster<W>::get_binary_bottom_cluster() const {
  return binary_cluster[1];
}

// Set unary cluster by idx
template <class W>
inline void Cluster<W>::set_unary_cluster(int idx, Cluster<W> *c) {
  unary_cluster[idx] = c;
}

// Get unary cluster by idx
template <class W>
inline Cluster<W>* Cluster<W>::get_unary_cluster(int idx) const {
  return unary_cluster[idx];
}

// Get maintained value
template <class W>
inline W Cluster<W>::get_val() const {
  return val;
}

// Set maintained value
template <class W>
inline void Cluster<W>::set_val(const W& rval) {
  val = rval;
}

// Get level
template <class W>
inline int Cluster<W>::get_level() const {
  return level;
}

// Get tid
template <class W>
inline uint64_t Cluster<W>::get_tid() const {
  return tid;
}

// Walking up the tree to set level and tid
template <class W>
void Cluster<W>::walk() {
  auto c = this;
  level = 0;
  tid = 0;
  while (c -> parent_cluster) {
    int pc = c -> parent_child_relation_by_pointer();
    c = c -> parent_cluster;
    ++level;
    ++tid;
    if (pc > 0) tid += c -> representative_cluster -> size;
    if (pc > 1 && c -> binary_cluster[0])  tid += c -> binary_cluster[0] -> size;
    if (pc > 2 && c -> binary_cluster[1])  tid += c -> binary_cluster[1] -> size;
    if (pc > 3 && c -> unary_cluster[0]) tid += c -> unary_cluster[0] -> size;
  }
}

// Set level and tid from parent, only for leaf clusters
template <class W>
inline void Cluster<W>::walk_from_parent_cluster() {
  level = parent_cluster -> level + 1;
  if (cluster_type == 0) {
    tid = parent_cluster -> tid + 1;
  } else {
    if (parent_child_relation_by_pointer() == 1)
      tid = parent_cluster -> tid + 2;
    else
      tid = parent_cluster -> tid + 3;
  }
}

// Print debug info for the rctree structure through dfs
template <class W>
void Cluster<W>::print() const {
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