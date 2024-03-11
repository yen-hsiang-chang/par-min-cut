#pragma once

#include <utility>

#include "parlay/sequence.h"
// They forgot to put pragma once.
#include "../pbbsbench/benchmarks/comparisonSort/ips4o/sort.h"

#define uint unsigned int

namespace utils {

template <class T, class BinPred>
void general_sort(parlay::sequence<T> &A, const BinPred& f) {
  ips4o::parallel::sort(A.begin(), A.end(), f);
}

template <class T>
T ternary_tree(const T& n, const parlay::sequence<T>& edge_ptr, 
               const parlay::sequence<std::tuple<T, T, T>>& edges,
               T& n_ternary, parlay::sequence<T>& edge_ptr_ternary, 
               parlay::sequence<std::tuple<T, T, T>>& edges_ternary) {
  
  auto deg = parlay::sequence<T>::from_function(
    n, [&](const T& i) {return edge_ptr[i + 1] - edge_ptr[i];});

  auto map_ternary = parlay::sequence<T>::from_function(
    n, [&](const T& i) {return deg[i] <= 3 ? 0 : deg[i] - 3;});

  n_ternary = parlay::scan_inplace(map_ternary) + n;
  map_ternary.push_back(n_ternary - n);
  
  edges_ternary.resize(2 * (n_ternary - 1));
  edge_ptr_ternary.resize(n_ternary);

  parlay::parallel_for(0, n, [&](const T& i) {
    if (deg[i] <= 3) {
      edge_ptr_ternary[i] = deg[i];
    } else {
      edge_ptr_ternary[i] = 3;
      parlay::parallel_for(n + map_ternary[i], n + map_ternary[i + 1], [&](const T& idx) {
        edge_ptr_ternary[idx] = 3;
      });
    }
  });

  edge_ptr_ternary.push_back(parlay::scan_inplace(edge_ptr_ternary));

  parlay::sequence<std::pair<T, T>> map_edge(n - 1);
  parlay::parallel_for(0, edges.size(), [&](const T& i) {
    auto [u, v, id] = edges[i];
    T offset = i - edge_ptr[u];
    T pos = deg[u] <= 3 ? u : (offset > 1 ? n + map_ternary[u] + offset - (offset + 1 == deg[u] ? 3 : 2) : u);
    if (u < v)
      map_edge[id].first = pos;
    else
      map_edge[id].second = pos;
  });

  parlay::parallel_for(0, n, [&](const T& i) {
    auto get_map_edge = [&](const std::tuple<T, T, T>& e) {
      auto [u, v, id] = e;
      auto res = map_edge[id];
      if(u >= v)  std::swap(res.first, res.second);
      return res;
    };
    if (deg[i] <= 3) {
      for (T idx = 0; idx < deg[i]; ++idx) {
        auto [mu, mv] = get_map_edge(edges[edge_ptr[i] + idx]);
        edges_ternary[edge_ptr_ternary[i] + idx] = std::make_tuple(mu, mv, std::get<2>(edges[edge_ptr[i] + idx]));
      }
    } else {
      {
        for (T idx = 0; idx < 2; ++idx) {
          auto [mu, mv] = get_map_edge(edges[edge_ptr[i] + idx]);
          edges_ternary[edge_ptr_ternary[i] + idx] = std::make_tuple(mu, mv, std::get<2>(edges[edge_ptr[i] + idx]));
        }
        edges_ternary[edge_ptr_ternary[i] + 2] = std::make_tuple(i, n + map_ternary[i], n - 1 + map_ternary[i]);
      }
      parlay::parallel_for(n + map_ternary[i], n + map_ternary[i + 1], [&](const T& idx) {
        auto [mu, mv] = get_map_edge(edges[edge_ptr[i] + idx - n - map_ternary[i] + 2]);
        edges_ternary[edge_ptr_ternary[idx]] = std::make_tuple(mu, mv, std::get<2>(edges[edge_ptr[i] + idx - n - map_ternary[i] + 2]));
        edges_ternary[edge_ptr_ternary[idx] + 1] = std::make_tuple(idx, idx - 1, idx - 1);
        edges_ternary[edge_ptr_ternary[idx] + 2] = std::make_tuple(idx, idx + 1, idx);
      });
      {
        edges_ternary[edge_ptr_ternary[n + map_ternary[i]] + 1] = std::make_tuple(n + map_ternary[i], i, n - 1 + map_ternary[i]);
        auto &e0 = edges_ternary[edge_ptr_ternary[n + map_ternary[i]]];
        auto &e1 = edges_ternary[edge_ptr_ternary[n + map_ternary[i]] + 1];
        if (std::get<1>(e0) > std::get<1>(e1))
          std::swap(e0, e1);
      }
      {
        auto [mu, mv] = get_map_edge(edges[edge_ptr[i + 1] - 1]);
        edges_ternary[edge_ptr_ternary[n + map_ternary[i + 1] - 1] + 2] = std::make_tuple(mu, mv, std::get<2>(edges[edge_ptr[i + 1] - 1]));
        auto &e1 = edges_ternary[edge_ptr_ternary[n + map_ternary[i + 1] - 1] + 1];
        auto &e2 = edges_ternary[edge_ptr_ternary[n + map_ternary[i + 1] - 1] + 2];
        if (std::get<1>(e1) > std::get<1>(e2))
          std::swap(e1, e2);
      }
    }
  });

  return parlay::find(deg, 1) - deg.begin();
}

template <class T>
void binary_tree(const T& root, const T& n_ternary, const parlay::sequence<T>& edge_ptr_ternary, 
                 const parlay::sequence<std::tuple<T, T, T>>& edges_ternary,
                 T& n_binary, parlay::sequence<T>& child_ptr_binary, 
                 parlay::sequence<std::pair<T, T>>& child_edges_binary,
                 parlay::sequence<std::pair<T, T>>& parent_binary) {
  
  auto deg_ternary = parlay::sequence<T>::from_function(
    n_ternary, [&](const T& i) {return edge_ptr_ternary[i + 1] - edge_ptr_ternary[i];});

  n_binary = n_ternary;

  child_ptr_binary = parlay::sequence<T>::from_function(
    n_binary, [&](const T& i) {return i == root ? deg_ternary[i] : deg_ternary[i] - 1;});

  child_ptr_binary.push_back(parlay::scan_inplace(child_ptr_binary));
  child_edges_binary.resize(n_binary - 1);
  parent_binary.resize(n_binary);
  
  parlay::sequence<T> val(edges_ternary.size()), nxtval(edges_ternary.size()), minval(edges_ternary.size());
  parlay::sequence<T> prv(edges_ternary.size()), nxtprv(edges_ternary.size());
  
  parlay::parallel_for(0, edges_ternary.size(), [&](const T& i) {
    auto [u, v, id] = edges_ternary[i];
    val[i] = 1;
    for(T j = edge_ptr_ternary[v]; j < edge_ptr_ternary[v + 1]; ++j) {
      if (std::get<1>(edges_ternary[j]) == u) {
        prv[j + 1 == edge_ptr_ternary[v + 1] ? edge_ptr_ternary[v] : j + 1] = i; 
        break;
      }
    }
  });
  prv[edge_ptr_ternary[root]] = edges_ternary.size();
  for (T step = 1; step < edges_ternary.size(); step <<= 1) {
    parlay::parallel_for(0, edges_ternary.size(), [&](const T& i) {
      if (prv[i] != edges_ternary.size()) {
        nxtprv[i] = prv[prv[i]];
        nxtval[i] = val[i] + val[prv[i]];
      } else {
        nxtprv[i] = prv[i];
        nxtval[i] = val[i];
      }
    });
    val.swap(nxtval);
    prv.swap(nxtprv);
  }
  parlay::parallel_for(0, n_ternary, [&](const T& i) {
    T res = edges_ternary.size();
    for(T j = edge_ptr_ternary[i]; j < edge_ptr_ternary[i + 1]; ++j) {
      res = std::min(res, val[j]);
    }
    minval[i] = res;
  });
  
  parlay::parallel_for(0, n_ternary, [&](const T& i) {
    T idx = child_ptr_binary[i];
    for(T j = edge_ptr_ternary[i]; j < edge_ptr_ternary[i + 1]; ++j) {
      auto [u, v, id] = edges_ternary[j];
      if (minval[u] < minval[v]) {
        child_edges_binary[idx++] = std::make_pair(v, id);
      } else {
        parent_binary[i] = std::make_pair(v, id);
      }
    }
  });
  parent_binary[root] = std::make_pair(n_binary, n_binary - 1);
}

} // namespace utils
