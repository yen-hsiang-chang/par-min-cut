#pragma once

#include <iostream>
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
void ternary_tree(const T& n, const T& m, const parlay::sequence<T>& edge_ptr, 
                  const parlay::sequence<std::tuple<T, T, T>>& edges, const parlay::sequence<T>& deg,
                  T& n_ternary, T& m_ternary, parlay::sequence<T>& edge_ptr_ternary, 
                  parlay::sequence<std::tuple<T, T, T>>& edges_ternary, parlay::sequence<T>& deg_ternary) {
  
  auto map_ternary = parlay::sequence<T>::from_function(
    n, [&](const T& i) {return deg[i] <= 3 ? 0 : deg[i] - 3;});
  n_ternary = parlay::scan_inplace(map_ternary) + n;
  map_ternary.push_back(n_ternary - n);
  m_ternary = 2 * (n_ternary - 1);
  
  edges_ternary.resize(m_ternary);
  deg_ternary.resize(n_ternary);

  parlay::parallel_for(0, n, [&](const T& i) {
    if (deg[i] <= 3) {
      deg_ternary[i] = deg[i];
    } else {
      deg_ternary[i] = 3;
      parlay::parallel_for(n + map_ternary[i], n + map_ternary[i + 1], [&](const T& idx) {
        deg_ternary[idx] = 3;
      });
    }
  });

  edge_ptr_ternary = deg_ternary;
  edge_ptr_ternary.push_back(parlay::scan_inplace(edge_ptr_ternary));

  parlay::sequence<std::pair<T, T>> map_edge(m / 2);
  parlay::parallel_for(0, m, [&](const T& i) {
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
        edges_ternary[edge_ptr_ternary[i] + 2] = std::make_tuple(i, n + map_ternary[i], m / 2 + map_ternary[i]);
      }
      parlay::parallel_for(n + map_ternary[i], n + map_ternary[i + 1], [&](const T& idx) {
        auto [mu, mv] = get_map_edge(edges[edge_ptr[i] + idx - n - map_ternary[i] + 2]);
        edges_ternary[edge_ptr_ternary[idx]] = std::make_tuple(mu, mv, std::get<2>(edges[edge_ptr[i] + idx - n - map_ternary[i] + 2]));
        edges_ternary[edge_ptr_ternary[idx] + 1] = std::make_tuple(idx, idx - 1, m / 2 + idx - n);
        edges_ternary[edge_ptr_ternary[idx] + 2] = std::make_tuple(idx, idx + 1, m / 2 + idx - n + 1);
      });
      {
        edges_ternary[edge_ptr_ternary[n + map_ternary[i]] + 1] = std::make_tuple(n + map_ternary[i], i, m / 2 + map_ternary[i]);
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
}

} // namespace utils
