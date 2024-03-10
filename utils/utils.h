#pragma once

#include "parlay/sequence.h"
// They forgot to put pragma once.
#include "../pbbsbench/benchmarks/comparisonSort/ips4o/sort.h"

#define uint unsigned int

namespace utils {

template <class T, class BinPred>
void general_sort(parlay::sequence<T> &A, const BinPred& f) {
  ips4o::parallel::sort(A.begin(), A.end(), f);
}

} // namespace utils
