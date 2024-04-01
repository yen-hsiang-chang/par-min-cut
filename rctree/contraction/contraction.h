#pragma once

#include <limits>
#include <utility>

#include "../cluster.h"
#include "../rctree.h"

struct Contraction_Index {
  size_t index[5];
  Contraction_Index(size_t i0 = 0, size_t i1 = 0, size_t i2 = 0, size_t i3 = 0, size_t i4 = 0) {
    index[0] = i0, index[1] = i1, index[2] = i2, index[3] = i3, index[4] = i4;
  }
  void set(const size_t& i) {index[i] = 1;}
  Contraction_Index operator+(const Contraction_Index& rhs) const {
    return Contraction_Index(index[0] + rhs[0], index[1] + rhs[1], index[2] + rhs[2], 
                             index[3] + rhs[3], index[4] + rhs[4]);
  }
  size_t operator[](const size_t& i) const {
    return index[i];
  }
};

template <class W>
struct Contraction_Type {
  W t, b; // sum of vertex weights of connected componenent from top/bottom boundaries, only count for current cluster 
  uint mt, mb; // max edge id from representative vertex to top/bottom boundaries
  bool joined;
  Contraction_Type(W t = 0, W b = 0, uint mt = 0, uint mb = 0, bool joined = false) :
                   t(t), b(b), mt(mt), mb(mb), joined(joined) {}
  Contraction_Type<W> operator+(const Contraction_Type<W>& rhs) const {
    return Contraction_Type<W>(t - rhs.t, b - rhs.b, 0, 0, joined || rhs.joined);
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
  // Constructor that first invokes base constructor then call build since build relies on virtual functions
  Contraction_RCTree(T n, const parlay::sequence<std::pair<T, T>>& edge_list,
                     parlay::random_generator& gen, 
                     const parlay::sequence<Contraction_Type<W>>& vertex_value,
                     const parlay::sequence<Contraction_Type<W>>& edge_value);
  // Virtual functions that are used for maintaining values in nullary/unary/binary clusters
  Contraction_Type<W> f_nullary(Cluster<T, Contraction_Type<W>> *c);
  Contraction_Type<W> f_unary(Cluster<T, Contraction_Type<W>> *c);
  Contraction_Type<W> f_binary(Cluster<T, Contraction_Type<W>> *c);
  // Base values for added vertex/edge when converting to ternary tree and rooted binary tree
  Contraction_Type<W> vertex_base();
  Contraction_Type<W> edge_base();

  // Find the maximum edge id between vertex u and v
  T query_max_path(T u, T v);
  // Sequential batch operations for testing
  W batch_operations_sequential(parlay::sequence<Contraction_MixOp<W>>& ops);
  W batch_operations(parlay::sequence<Contraction_MixOp<W>>& ops);

private:
  // Operation 0: Subtract vertex weight
  void subtract_vertex_weight(T u, W w);
  // Operation 1: Join edge
  void join_edge(T e);
  // Operation 2: Query vertex sum
  W query_vertex_sum(T u);
  // Extract maximum edge id from query to top/bottom boundaries
  void extract_max_path(Cluster<T, Contraction_Type<W>> *c, std::pair<T, T>& mx);
};

// Constructor that first invokes base constructor then call build since build relies on virtual functions
template <class T, class W>
Contraction_RCTree<T, W>::Contraction_RCTree(T n, const parlay::sequence<std::pair<T, T>>& edge_list,
                                             parlay::random_generator& gen, 
                                             const parlay::sequence<Contraction_Type<W>>& vertex_value,
                                             const parlay::sequence<Contraction_Type<W>>& edge_value)
                                             : RCTree<T, Contraction_Type<W>>(n, edge_list) {
  this -> build(gen, vertex_value, edge_value);
}

// Virtual function that is used for maintaining values in nullary clusters
template <class T, class W>
Contraction_Type<W> Contraction_RCTree<T, W>::f_nullary(Cluster<T, Contraction_Type<W>> *c) {
  return Contraction_Type<W>(0, 0, 0, 0, true);
}

// Virtual function that is used for maintaining values in unary clusters
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

// Virtual function that is used for maintaining values in binary clusters
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

// Base value for added vertex when converting to ternary tree and rooted binary tree
template <class T, class W>
inline Contraction_Type<W> Contraction_RCTree<T, W>::vertex_base() {
  return Contraction_Type<W>(0, 0, 0, 0, true);
}

// Base value for added edge when converting to ternary tree and rooted binary tree
template <class T, class W>
inline Contraction_Type<W> Contraction_RCTree<T, W>::edge_base() {
  return Contraction_Type<W>(0, 0, 0, 0, true);
}

// Operation 0: Subtract vertex weight
template <class T, class W>
void Contraction_RCTree<T, W>::subtract_vertex_weight(T u, W w) {
  auto val = RCTree<T, Contraction_Type<W>>::vertex_clusters[u].get_val();
  val.t -= w;
  RCTree<T, Contraction_Type<W>>::vertex_clusters[u].set_val(val);
  RCTree<T, Contraction_Type<W>>::reevaluate(&RCTree<T, Contraction_Type<W>>::vertex_clusters[u]);
}

// Operation 1: Join edge
template <class T, class W>
void Contraction_RCTree<T, W>::join_edge(T e) {
  auto val = RCTree<T, Contraction_Type<W>>::edge_clusters[e].get_val();
  val.joined = true;
  RCTree<T, Contraction_Type<W>>::edge_clusters[e].set_val(val);
  RCTree<T, Contraction_Type<W>>::reevaluate(&RCTree<T, Contraction_Type<W>>::edge_clusters[e]);
}

// Operation 2: Query vertex sum
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

// Extract maximum edge id from query to top/bottom boundaries
template <class T, class W>
void Contraction_RCTree<T, W>::extract_max_path(Cluster<T, Contraction_Type<W>> *c, std::pair<T, T>& mx) {
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

// Find the maximum edge id between vertex u and v
template <class T, class W>
T Contraction_RCTree<T, W>::query_max_path(T u, T v) {
  auto cu = &RCTree<T, Contraction_Type<W>>::vertex_clusters[u];
  auto cv = &RCTree<T, Contraction_Type<W>>::vertex_clusters[v];
  T du = cu -> get_level(), dv = cv -> get_level();
  if (du < dv)
    std::swap(cu, cv), std::swap(du, dv);
  std::pair<T, T> mu = {0, 0}, mv = {0, 0}; 
  while(du > dv)  --du, extract_max_path(cu, mu), cu = cu -> get_parent_cluster();

  while(cu -> get_parent_cluster() != cv -> get_parent_cluster()) {
    extract_max_path(cu, mu), extract_max_path(cv, mv);
    cu = cu -> get_parent_cluster(), cv = cv -> get_parent_cluster();
  }
  std::pair<T, T> bu, bv;
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

template <class T, class W>
W Contraction_RCTree<T, W>::batch_operations_sequential(parlay::sequence<Contraction_MixOp<W>>& ops) {
  utils::general_sort(ops, [&](const Contraction_MixOp<W>& lhs, const Contraction_MixOp<W>& rhs) {
    return lhs.timestamp < rhs.timestamp;
  });
  W answer = std::numeric_limits<W>::max();
  for (auto q : ops) {
    if (q.type == 0) {
      subtract_vertex_weight(q.u, q.w);
    } else if(q.type == 1) {
      join_edge(q.u);
    } else {
      answer = std::min(answer, query_vertex_sum(q.u));
    }
  }
  return answer;
}

template <class T, class W>
W Contraction_RCTree<T, W>::batch_operations(parlay::sequence<Contraction_MixOp<W>>& ops) {
  utils::general_sort(ops, [&](const Contraction_MixOp<W>& lhs, const Contraction_MixOp<W>& rhs) {
    T ltid = lhs.type == 1 ? RCTree<T, Contraction_Type<W>>::edge_clusters[lhs.u].get_tid() 
                           : RCTree<T, Contraction_Type<W>>::vertex_clusters[lhs.u].get_tid();
    T rtid = rhs.type == 1 ? RCTree<T, Contraction_Type<W>>::edge_clusters[rhs.u].get_tid()
                           : RCTree<T, Contraction_Type<W>>::vertex_clusters[rhs.u].get_tid();
    return ltid == rtid ? lhs.timestamp < rhs.timestamp : ltid < rtid;
  });

  parlay::sequence<std::pair<size_t, Cluster<T, Contraction_Type<W>>>> current_level(ops.size());
  parlay::sequence<std::pair<Contraction_Type<W>, bool>> segmented_contraction(ops.size());
  parlay::parallel_for(0, ops.size(), [&](const size_t& i) {
    current_level[i] = std::make_pair(i, ops[i].type == 1 ? RCTree<T, Contraction_Type<W>>::edge_clusters[ops[i].u]
                                                          : RCTree<T, Contraction_Type<W>>::vertex_clusters[ops[i].u]);
  });
  parlay::parallel_for(0, ops.size(), [&](const size_t& i) {
    if (i == 0 || current_level[i].second.get_tid() != current_level[i - 1].second.get_tid()) {
      Contraction_Type<W> rhs(ops[i].type == 0 ? ops[i].w : 0, 0, 0, 0, ops[i].type == 1 ? true : false);
      segmented_contraction[i].first = current_level[i].second.get_val() + rhs;
      segmented_contraction[i].second = true;
    } else {
      Contraction_Type<W> rhs(ops[i].type == 0 ? ops[i].w : 0, 0, 0, 0, ops[i].type == 1 ? true : false);
      segmented_contraction[i].first = rhs;
      segmented_contraction[i].second = false;
    }
  });

  auto segmented_contraction_scan_f = [](const std::pair<Contraction_Type<W>, bool>& lhs,
                                         const std::pair<Contraction_Type<W>, bool>& rhs) {
    return std::make_pair(rhs.second ? rhs.first : lhs.first + rhs.first, lhs.second || rhs.second);
  };
  parlay::scan_inclusive_inplace(segmented_contraction, parlay::make_monoid(segmented_contraction_scan_f, 
                                 std::make_pair(Contraction_Type<W>(), false)));
  parlay::parallel_for(0, ops.size(), [&](const size_t& i) {
    current_level[i].second.set_val(segmented_contraction[i].first);
  });
  auto max_level_seq = parlay::delayed_tabulate(ops.size(), [&](const size_t& i) {
    return current_level[i].second.get_level();
  });
  int max_level = parlay::reduce(max_level_seq, parlay::maximum<int>());
  
  parlay::sequence<std::pair<size_t, bool>> segmented_level(ops.size());
  auto segmented_level_scan_f = [](const std::pair<size_t, bool>& lhs,
                                   const std::pair<size_t, bool>& rhs) {
    return std::make_pair(rhs.second ? rhs.first : lhs.first + rhs.first, lhs.second || rhs.second);
  };

  for (int level = max_level; level >= 1; level--) {
    parlay::parallel_for(0, ops.size(), [&](const size_t& i) {
      segmented_level[i] = std::make_pair(
        current_level[i].second.get_level() == level ? 1 : 0,
        i == 0 || current_level[i].second.get_tid() != current_level[i - 1].second.get_tid() ? true : false
      );
    });
    parlay::scan_inclusive_inplace(segmented_level, parlay::make_monoid(segmented_level_scan_f, 
                                   std::make_pair(size_t(0), false)));
    parlay::parallel_for(0, ops.size(), [&](const size_t& i) {
      segmented_level[i].second = false;
      if (segmented_level[i].first != 0) {
        if (i + 1 == ops.size() || (current_level[i].second.get_parent_cluster() != current_level[i + 1].second.get_parent_cluster())) {
          segmented_level[i].second = true;
        }
      }
    });
    parlay::parallel_for(0, ops.size(), [&](const size_t& i) {
      if (segmented_level[i].second) {
        size_t idx = i;
        auto merged = current_level.subseq(idx + 1 - segmented_level[idx].first, idx + 1);
        size_t offset[5] = {0, 0, 0, 0, 0};
        offset[current_level[idx].second.parent_child_relation_by_tid()] = idx + 1 - segmented_level[idx].first;
        while (idx >= segmented_level[idx].first) {
          idx -= segmented_level[idx].first;
          if (segmented_level[idx].first == 0 || current_level[idx].second.get_parent_cluster() != current_level[idx + 1].second.get_parent_cluster())
            break;
          merged = parlay::merge(parlay::make_slice(merged), parlay::make_slice(current_level.begin() + idx + 1 - segmented_level[idx].first, current_level.begin() + idx + 1), 
                                 [&](const std::pair<size_t, Cluster<T, Contraction_Type<W>>>& lhs,
                                     const std::pair<size_t, Cluster<T, Contraction_Type<W>>>& rhs) {
            return ops[lhs.first].timestamp < ops[rhs.first].timestamp;
          });
          offset[current_level[idx].second.parent_child_relation_by_tid()] = idx + 1 - segmented_level[idx].first;
        }
        parlay::sequence<Contraction_Index> cindex(merged.size());
        parlay::parallel_for(0, merged.size(), [&](const size_t& j) {
          cindex[j].set(merged[j].second.parent_child_relation_by_tid());
        });
        parlay::scan_inclusive_inplace(cindex);
        parlay::parallel_for(0, merged.size(), [&](const size_t& j) {
          auto &c = merged[j].second;
          c = *c.get_parent_cluster();
          if (cindex[j][0]) {
            c.set_representative_cluster(&current_level[offset[0] + cindex[j][0] - 1].second);
          }
          if (cindex[j][1]) {
            c.set_binary_top_cluster(&current_level[offset[1] + cindex[j][1] - 1].second);
          }
          if (cindex[j][2]) {
            c.set_binary_bottom_cluster(&current_level[offset[2] + cindex[j][2] - 1].second);
          }
          if (cindex[j][3]) {
            c.set_unary_cluster(0, &current_level[offset[3] + cindex[j][3] - 1].second);
          }
          if (cindex[j][4]) {
            c.set_unary_cluster(1, &current_level[offset[4] + cindex[j][4] - 1].second);
          }
          switch(c.get_cluster_type()) {
            case 0: c.set_val(f_nullary(&c)); break;
            case 1: c.set_val(f_unary(&c)); break;
            case 2: c.set_val(f_binary(&c)); break;
            default: assert(false);
          }
          auto &q = ops[merged[j].first];
          if (q.type == 2 && c.get_representative() == q.u) {
            if (c.get_binary_top_cluster() && c.get_binary_top_cluster() -> get_val().joined) {
              q.u = c.get_top_boundary();
            } else if (c.get_binary_bottom_cluster() && c.get_binary_bottom_cluster() -> get_val().joined) {
              q.u = c.get_bottom_boundary();
            } else {
              W qsum = c.get_representative_cluster() -> get_val().t;
              for (int k = 0; k < 2; ++k)
                if (c.get_unary_cluster(k))
                  qsum += c.get_unary_cluster(k) -> get_val().t;
              if (c.get_binary_top_cluster())
                qsum += c.get_binary_top_cluster() -> get_val().b;
              if (c.get_binary_bottom_cluster())
                qsum += c.get_binary_bottom_cluster() -> get_val().t;
              q.w = qsum;
              q.type = 3;
            }
          }
        });
        parlay::copy(merged, current_level.cut(i + 1 - merged.size(), i + 1));
      }
    });
  }
  auto min_cut_seq = parlay::delayed_tabulate(ops.size(), [&](const size_t& i) {
    return ops[i].type == 3 ? ops[i].w : std::numeric_limits<W>::max();
  });
  return parlay::reduce(min_cut_seq, parlay::minimum<int>());;
}