#pragma once

#include <limits>
#include <utility>

#include "../cluster.h"
#include "../rctree.h"
#include "../utils/utils.h"

struct Contraction_Index {
  uint64_t index[5];
  Contraction_Index(uint64_t i0 = 0, uint64_t i1 = 0, uint64_t i2 = 0, uint64_t i3 = 0, uint64_t i4 = 0) {
    index[0] = i0, index[1] = i1, index[2] = i2, index[3] = i3, index[4] = i4;
  }
  void set(int i) {index[i] = 1;}
  Contraction_Index operator+(const Contraction_Index& rhs) const {
    return Contraction_Index(index[0] + rhs[0], index[1] + rhs[1], index[2] + rhs[2], 
                             index[3] + rhs[3], index[4] + rhs[4]);
  }
  uint64_t operator[](int i) const {
    return index[i];
  }
};

struct Contraction_P {
  uintV mt, mb; // max edge id from representative vertex to top/bottom boundaries
  Contraction_P(uintV mt = 0, uintV mb = 0) : mt(mt), mb(mb) {}
};

template <class W>
struct Contraction_CC {
  W t, b;
  bool joined;
  Contraction_CC(W t = 0, W b = 0, bool joined = false) :
                 t(t), b(b), joined(joined) {}
  Contraction_CC<W> operator+(const Contraction_CC<W>& rhs) const {
    return Contraction_CC<W>(t + rhs.t, b + rhs.b, joined || rhs.joined);
  }
};

template <class W>
using Contraction_Type = std::pair<Contraction_CC<W>, Contraction_P>;

template <class W>
struct Contraction_MixOp {
  uint64_t timestamp, tid;
  int type;
  uintV u;
  W w;
  Contraction_MixOp(uint64_t timestamp = 0, int type = 0, uintV u = 0, W w = 0) :
                    timestamp(timestamp), type(type), u(u), w(w) {}
};

template <class W>
class Contraction_RCTree : public RCTree<Contraction_Type<W>> {
public:
  // Constructor that first invokes base constructor then call build since build relies on virtual functions
  Contraction_RCTree(uintV n, const parlay::sequence<std::pair<uintV, uintV>>& edge_list,
                     parlay::random_generator& gen, 
                     const parlay::sequence<Contraction_Type<W>>& vertex_value,
                     const parlay::sequence<Contraction_Type<W>>& edge_value);
  // Virtual functions that are used for maintaining values in nullary/unary/binary clusters
  Contraction_Type<W> f_nullary(Cluster<Contraction_Type<W>> *c);
  Contraction_Type<W> f_unary(Cluster<Contraction_Type<W>> *c);
  Contraction_Type<W> f_binary(Cluster<Contraction_Type<W>> *c);
  // Base values for added vertex/edge when converting to ternary tree and rooted binary tree
  Contraction_Type<W> vertex_base();
  Contraction_Type<W> edge_base();

  // Find the maximum edge id between vertex u and v
  uintV query_max_path(uintV u, uintV v);
  // Sequential batch operations for testing
  W batch_operations_sequential(parlay::sequence<Contraction_MixOp<W>>& ops);
  W batch_operations(parlay::sequence<Contraction_MixOp<W>>& ops);

private:
  // Operation 0: Subtract vertex weight
  void subtract_vertex_weight(uintV u, W w);
  // Operation 1: Join edge
  void join_edge(uintV e);
  // Operation 2: Query vertex sum
  W query_vertex_sum(uintV u);
  // Extract maximum edge id from query to top/bottom boundaries
  void extract_max_path(Cluster<Contraction_Type<W>> *c, std::pair<uintV, uintV>& mx);
    Contraction_CC<W> f_nullary_CC() {
    return Contraction_CC<W>(0, 0, true);
  }
  Contraction_CC<W> f_unary_CC(W s, Contraction_CC<W>& binary_top_val) {
    if (!binary_top_val.joined)
      return Contraction_CC<W>(binary_top_val.t, 0, true);
    return Contraction_CC<W>(s + binary_top_val.t, 0, true);
  }
  Contraction_CC<W> f_binary_CC(W s, Contraction_CC<W>& binary_top_val, Contraction_CC<W>& binary_bottom_val) {
    if (binary_top_val.joined) {
      if (binary_bottom_val.joined)
        return Contraction_CC<W>(s + binary_top_val.t + binary_bottom_val.t, 0, true);
      return Contraction_CC<W>(s + binary_top_val.t + binary_bottom_val.t, binary_bottom_val.b, false);
    }
    if (binary_bottom_val.joined)
      return Contraction_CC<W>(binary_top_val.t, s + binary_top_val.b + binary_bottom_val.t, false);
    return Contraction_CC<W>(binary_top_val.t, binary_bottom_val.b, false);
  }
};

// Constructor that first invokes base constructor then call build since build relies on virtual functions
template <class W>
Contraction_RCTree<W>::Contraction_RCTree(uintV n, const parlay::sequence<std::pair<uintV, uintV>>& edge_list,
                                          parlay::random_generator& gen, 
                                          const parlay::sequence<Contraction_Type<W>>& vertex_value,
                                          const parlay::sequence<Contraction_Type<W>>& edge_value)
                                          : RCTree<Contraction_Type<W>>(n, edge_list) {
  this -> build(gen, vertex_value, edge_value);
}

// Virtual function that is used for maintaining values in nullary clusters
template <class W>
Contraction_Type<W> Contraction_RCTree<W>::f_nullary(Cluster<Contraction_Type<W>> *c) {
  return std::make_pair(Contraction_CC<W>(0, 0, true), Contraction_P(0, 0));
}

// Virtual function that is used for maintaining values in unary clusters
template <class W>
Contraction_Type<W> Contraction_RCTree<W>::f_unary(Cluster<Contraction_Type<W>> *c) {
  Contraction_P cp(std::max(c -> get_binary_top_cluster() -> get_val().second.mt, c -> get_binary_top_cluster() -> get_val().second.mb), 0);
  if (!c -> get_binary_top_cluster() -> get_val().first.joined)  
    return std::make_pair(Contraction_CC<W>(c -> get_binary_top_cluster() -> get_val().first.t, 0, true), cp);
  W s = c -> get_representative_cluster() -> get_val().first.t;
  for (int i = 0; i < 2; ++i)
    if (c -> get_unary_cluster(i))
      s += c -> get_unary_cluster(i) -> get_val().first.t;
  return std::make_pair(Contraction_CC<W>(s + c -> get_binary_top_cluster() -> get_val().first.t, 0, true), cp);
}

// Virtual function that is used for maintaining values in binary clusters
template <class W>
Contraction_Type<W> Contraction_RCTree<W>::f_binary(Cluster<Contraction_Type<W>> *c) {
  Contraction_P cp(std::max(c -> get_binary_top_cluster() -> get_val().second.mt, c -> get_binary_top_cluster() -> get_val().second.mb), 
                   std::max(c -> get_binary_bottom_cluster() -> get_val().second.mt, c -> get_binary_bottom_cluster() -> get_val().second.mb));
  W s = c -> get_representative_cluster() -> get_val().first.t;
  for (int i = 0; i < 2; ++i)
    if (c -> get_unary_cluster(i))
      s += c -> get_unary_cluster(i) -> get_val().first.t;
  if (c -> get_binary_top_cluster() -> get_val().first.joined) {
    if (c -> get_binary_bottom_cluster() -> get_val().first.joined)
      return std::make_pair(Contraction_CC<W>(s + c -> get_binary_top_cluster() -> get_val().first.t + c -> get_binary_bottom_cluster() -> get_val().first.t, 
                                              0, true), cp);
    return  std::make_pair(Contraction_CC<W>(c -> get_binary_top_cluster() -> get_val().first.t + c -> get_binary_bottom_cluster() -> get_val().first.t + s, 
                                             c -> get_binary_bottom_cluster() -> get_val().first.b, false), cp);
  }
  if (c -> get_binary_bottom_cluster() -> get_val().first.joined)
    return std::make_pair(Contraction_CC<W>(c -> get_binary_top_cluster() -> get_val().first.t, 
                                            c -> get_binary_top_cluster() -> get_val().first.b + c -> get_binary_bottom_cluster() -> get_val().first.t + s, false), cp);;
  return std::make_pair(Contraction_CC<W>(c -> get_binary_top_cluster() -> get_val().first.t, 
                                          c -> get_binary_bottom_cluster() -> get_val().first.b, false), cp);
}

// Base value for added vertex when converting to ternary tree and rooted binary tree
template <class W>
inline Contraction_Type<W> Contraction_RCTree<W>::vertex_base() {
  return std::make_pair(Contraction_CC<W>(0, 0, true), Contraction_P(0, 0));
}

// Base value for added edge when converting to ternary tree and rooted binary tree
template <class W>
inline Contraction_Type<W> Contraction_RCTree<W>::edge_base() {
  return std::make_pair(Contraction_CC<W>(0, 0, true), Contraction_P(0, 0));
}

// Operation 0: Subtract vertex weight
template <class W>
void Contraction_RCTree<W>::subtract_vertex_weight(uintV u, W w) {
  auto val = RCTree<Contraction_Type<W>>::vertex_clusters[u].get_val();
  val.first.t -= w;
  RCTree<Contraction_Type<W>>::vertex_clusters[u].set_val(val);
  RCTree<Contraction_Type<W>>::reevaluate(&RCTree<Contraction_Type<W>>::vertex_clusters[u]);
}

// Operation 1: Join edge
template <class W>
void Contraction_RCTree<W>::join_edge(uintV e) {
  auto val = RCTree<Contraction_Type<W>>::edge_clusters[e].get_val();
  val.first.joined = true;
  RCTree<Contraction_Type<W>>::edge_clusters[e].set_val(val);
  RCTree<Contraction_Type<W>>::reevaluate(&RCTree<Contraction_Type<W>>::edge_clusters[e]);
}

// Operation 2: Query vertex sum
template <class W>
W Contraction_RCTree<W>::query_vertex_sum(uintV u) {
  auto c = RCTree<Contraction_Type<W>>::vertex_clusters[u].get_parent_cluster();
  while (true) {
    if (c -> get_binary_top_cluster() && c -> get_binary_top_cluster() -> get_val().first.joined) {
      c = RCTree<Contraction_Type<W>>::vertex_clusters[c -> get_top_boundary()].get_parent_cluster();
      continue;
    }
    if (c -> get_binary_bottom_cluster() && c -> get_binary_bottom_cluster() -> get_val().first.joined) {
      c = RCTree<Contraction_Type<W>>::vertex_clusters[c -> get_bottom_boundary()].get_parent_cluster();
      continue;
    }
    break;
  }
  W s = c -> get_representative_cluster() -> get_val().first.t;
  for (int i = 0; i < 2; ++i)
    if (c -> get_unary_cluster(i))
      s += c -> get_unary_cluster(i) -> get_val().first.t;
  if (c -> get_binary_top_cluster())
    s += c -> get_binary_top_cluster() -> get_val().first.b;
  if (c -> get_binary_bottom_cluster())
    s += c -> get_binary_bottom_cluster() -> get_val().first.t;
  return s;
}

// Extract maximum edge id from query to top/bottom boundaries
template <class W>
void Contraction_RCTree<W>::extract_max_path(Cluster<Contraction_Type<W>> *c, std::pair<uintV, uintV>& mx) {
  if (c -> get_parent_cluster() -> get_cluster_type() == 1) {
    if (c -> get_cluster_type() == 2)
      mx = {mx.first, 0};
    else
      mx = {std::max(mx.first, c -> get_parent_cluster() -> get_val().second.mt), 0};
  } else {
    if (c -> get_cluster_type() == 2) {
      if (c -> get_parent_cluster() -> get_binary_top_cluster() == c)
        mx = {mx.first, std::max(mx.second, c -> get_parent_cluster() -> get_val().second.mb)};
      else
        mx = {std::max(mx.first, c -> get_parent_cluster() -> get_val().second.mt), mx.second};
    } else {
      mx = {std::max(mx.first, c -> get_parent_cluster() -> get_val().second.mt), std::max(mx.first, c -> get_parent_cluster() -> get_val().second.mb)};
    }
  }
}

// Find the maximum edge id between vertex u and v
template <class W>
uintV Contraction_RCTree<W>::query_max_path(uintV u, uintV v) {
  auto cu = &RCTree<Contraction_Type<W>>::vertex_clusters[u];
  auto cv = &RCTree<Contraction_Type<W>>::vertex_clusters[v];
  int du = cu -> get_level(), dv = cv -> get_level();
  if (du < dv)
    std::swap(cu, cv), std::swap(du, dv);
  std::pair<uintV, uintV> mu = {0, 0}, mv = {0, 0}; 
  while(du > dv)  --du, extract_max_path(cu, mu), cu = cu -> get_parent_cluster();

  while(cu -> get_parent_cluster() != cv -> get_parent_cluster()) {
    extract_max_path(cu, mu), extract_max_path(cv, mv);
    cu = cu -> get_parent_cluster(), cv = cv -> get_parent_cluster();
  }
  std::pair<uintV, uintV> bu, bv;
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

template <class W>
W Contraction_RCTree<W>::batch_operations_sequential(parlay::sequence<Contraction_MixOp<W>>& ops) {
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

template <class W>
W Contraction_RCTree<W>::batch_operations(parlay::sequence<Contraction_MixOp<W>>& ops) {
  {
    auto vertex_tid = parlay::sequence<uint64_t>::from_function(RCTree<Contraction_Type<W>>::n, [&](uintV i) {
      return RCTree<Contraction_Type<W>>::vertex_clusters[i].get_tid();
    });
    auto edge_tid = parlay::sequence<uint64_t>::from_function(RCTree<Contraction_Type<W>>::n - 1, [&](uintV i) {
      return RCTree<Contraction_Type<W>>::edge_clusters[i].get_tid();
    });
    parlay::parallel_for(0, ops.size(), [&](uint64_t i) {
      ops[i].tid = ops[i].type == 1 ? edge_tid[ops[i].u] : vertex_tid[ops[i].u];
    });
  }
  
  utils::general_sort(ops, [&](const Contraction_MixOp<W>& lhs, const Contraction_MixOp<W>& rhs) {
    return lhs.tid == rhs.tid ? lhs.timestamp < rhs.timestamp : lhs.tid < rhs.tid;
  });

  parlay::sequence<std::tuple<uint64_t, Cluster<Contraction_Type<W>>*, Contraction_CC<W>>> current_level(ops.size());
  parlay::sequence<std::pair<Contraction_CC<W>, bool>> segmented_contraction(ops.size());
  parlay::parallel_for(0, ops.size(), [&](uint64_t i) {
    std::get<0>(current_level[i]) = i;
    std::get<1>(current_level[i]) = ops[i].type == 1 ? &RCTree<Contraction_Type<W>>::edge_clusters[ops[i].u]
                                                     : &RCTree<Contraction_Type<W>>::vertex_clusters[ops[i].u];
  });
  
  parlay::parallel_for(0, ops.size(), [&](uint64_t i) {
    if (i == 0 || ops[i].tid != ops[i - 1].tid) {
      Contraction_CC<W> rhs(ops[i].type == 0 ? -ops[i].w : 0, 0, ops[i].type == 1 ? true : false);
      segmented_contraction[i].first = std::get<1>(current_level[i]) -> get_val().first + rhs;
      segmented_contraction[i].second = true;
    } else {
      Contraction_CC<W> rhs(ops[i].type == 0 ? -ops[i].w : 0, 0, ops[i].type == 1 ? true : false);
      segmented_contraction[i].first = rhs;
      segmented_contraction[i].second = false;
    }
  });

  auto segmented_contraction_scan_f = [](const std::pair<Contraction_CC<W>, bool>& lhs,
                                         const std::pair<Contraction_CC<W>, bool>& rhs) {
    return std::make_pair(rhs.second ? rhs.first : lhs.first + rhs.first, lhs.second || rhs.second);
  };
  parlay::scan_inclusive_inplace(segmented_contraction, parlay::make_monoid(segmented_contraction_scan_f, 
                                 std::make_pair(Contraction_CC<W>(), false)));

  parlay::parallel_for(0, ops.size(), [&](uint64_t i) {
    std::get<2>(current_level[i]) = segmented_contraction[i].first;
  });

  auto level_seq = parlay::sequence<int>::from_function(RCTree<Contraction_Type<W>>::n_binary, [&](uintV i) {
    return RCTree<Contraction_Type<W>>::vertex_clusters[i].get_level();
  });
  int max_level = parlay::reduce(level_seq, parlay::maximum<int>());
  parlay::sequence<uint64_t> pos_seq(RCTree<Contraction_Type<W>>::n_binary, ops.size());
  
  parlay::sequence<std::pair<uint64_t, bool>> segmented_level(ops.size());
  auto segmented_level_scan_f = [](const std::pair<uint64_t, bool>& lhs,
                                   const std::pair<uint64_t, bool>& rhs) {
    return std::make_pair(rhs.second ? rhs.first : lhs.first + rhs.first, lhs.second || rhs.second);
  };

  parlay::parallel_for(0, ops.size(), [&](uint64_t i) {
    segmented_level[i] = std::make_pair(1, i + 1 == ops.size() || ops[i].tid != ops[i + 1].tid);
  });
  parlay::scan_inclusive_inplace(parlay::make_slice(segmented_level.rbegin(), segmented_level.rend()), parlay::make_monoid(segmented_level_scan_f, 
                                  std::make_pair(uint64_t(0), false)));
  
  parlay::parallel_for(0, ops.size(), [&](uint64_t i) {
    if (i == 0 || (std::get<1>(current_level[i]) -> get_parent_cluster() -> get_tid() > ops[i - 1].tid)) {
      pos_seq[std::get<1>(current_level[i]) -> get_parent_cluster() -> get_representative()] = i;
      segmented_level[i].second = false;
    }
    else {
      segmented_level[i].second = true;
    }
  });

  for (int level = max_level; level >= 1; --level) {
    parlay::parallel_for(0, level_seq.size(), [&](uintV u) {
      if (level_seq[u] == level && pos_seq[u] != ops.size()) {
        uint64_t idx = pos_seq[u];
        auto merged = current_level.subseq(idx, idx + segmented_level[idx].first);
        uint64_t offset[5] = {0, 0, 0, 0, 0};
        offset[std::get<1>(current_level[idx]) -> parent_child_relation_by_tid()] = idx;
        while (idx + segmented_level[idx].first < ops.size()) {
          idx += segmented_level[idx].first;
          if (!segmented_level[idx].second || std::get<1>(current_level[idx]) -> get_parent_cluster() != std::get<1>(current_level[pos_seq[u]]) -> get_parent_cluster())
            break;
          merged = parlay::merge(parlay::make_slice(merged), current_level.cut(idx, idx + segmented_level[idx].first),
                                 [&](const std::tuple<uint64_t, Cluster<Contraction_Type<W>>*, Contraction_CC<W>>& lhs,
                                     const std::tuple<uint64_t, Cluster<Contraction_Type<W>>*, Contraction_CC<W>>& rhs) {
            return ops[std::get<0>(lhs)].timestamp < ops[std::get<0>(rhs)].timestamp;
          });
          offset[std::get<1>(current_level[idx]) -> parent_child_relation_by_tid()] = idx;
        }
        parlay::sequence<Contraction_Index> cindex(merged.size());
        parlay::parallel_for(0, merged.size(), [&](uint64_t j) {
          cindex[j].set(std::get<1>(merged[j]) -> parent_child_relation_by_tid());
        });
        parlay::scan_inclusive_inplace(cindex);
        parlay::parallel_for(0, merged.size(), [&](uint64_t j) {
          auto &c = std::get<1>(merged[j]);
          c = c -> get_parent_cluster();
          Contraction_CC<W> binary_top_val, binary_bottom_val;
          W s = 0;
          if (cindex[j][0]) {
            s += std::get<2>(current_level[offset[0] + cindex[j][0] - 1]).t;
          } else {
            s += c -> get_representative_cluster() -> get_val().first.t;
          }
          if (cindex[j][1]) {
            binary_top_val = std::get<2>(current_level[offset[1] + cindex[j][1] - 1]);
          } else {
            if (c -> get_cluster_type() > 0) {   
              binary_top_val = c -> get_binary_top_cluster() -> get_val().first;
            }
          }
          if (cindex[j][2]) {
            binary_bottom_val = std::get<2>(current_level[offset[2] + cindex[j][2] - 1]);
          } else {
            if (c -> get_cluster_type() > 1) {
              binary_bottom_val = c -> get_binary_bottom_cluster() -> get_val().first;
            }
          }
          if (cindex[j][3]) {
            s += std::get<2>(current_level[offset[3] + cindex[j][3] - 1]).t;
          } else {
            if (c -> get_unary_cluster(0)) {
              s += c -> get_unary_cluster(0) -> get_val().first.t;
            }
          }
          if (cindex[j][4]) {
            s += std::get<2>(current_level[offset[4] + cindex[j][4] - 1]).t;
          } else {
            if (c -> get_unary_cluster(1)) {
              s += c -> get_unary_cluster(1) -> get_val().first.t;
            }
          }
          switch(c -> get_cluster_type()) {
            case 0: std::get<2>(merged[j]) = f_nullary_CC(); break;
            case 1: std::get<2>(merged[j]) = f_unary_CC(s, binary_top_val); break;
            case 2: std::get<2>(merged[j]) = f_binary_CC(s, binary_top_val, binary_bottom_val); break;
            default: assert(false);
          }
          auto &q = ops[std::get<0>(merged[j])];
          if (q.type == 2 && c -> get_representative() == q.u) {
            if (binary_top_val.joined) {
              q.u = c -> get_top_boundary();
            } else if (binary_bottom_val.joined) {
              q.u = c -> get_bottom_boundary();
            } else {
              q.w = s + binary_top_val.b + binary_bottom_val.t;
              q.type = 3;
            }
          }
        });
        parlay::copy(merged, current_level.cut(pos_seq[u], pos_seq[u] + merged.size()));
        segmented_level[pos_seq[u]].first = merged.size();
      }
    });
    parlay::parallel_for(0, level_seq.size(), [&](uintV u) {
      if (level_seq[u] == level && pos_seq[u] != ops.size()) {
        auto i = pos_seq[u];
        if (i == 0 || std::get<1>(current_level[i]) -> get_parent_cluster() != std::get<1>(current_level[i - 1]) -> get_parent_cluster()) {
          if (std::get<1>(current_level[i]) -> get_parent_cluster())
            pos_seq[std::get<1>(current_level[i]) -> get_parent_cluster() -> get_representative()] = i;
        }
        else {
          segmented_level[i].second = true;
        }
      }
    });
  }
  auto min_cut_seq = parlay::delayed_tabulate(ops.size(), [&](uint64_t i) {
    return ops[i].type == 3 ? ops[i].w : std::numeric_limits<W>::max();
  });
  return parlay::reduce(min_cut_seq, parlay::minimum<W>());;
}