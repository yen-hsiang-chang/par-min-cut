# Parallel Randomized Minimum Cuts and Tree Contractions

## Description

This is a personal project trying to implement the first half of the parallel randomized minimum cut algorithm proposed in [Parallel Minimum Cuts in $O(m\log^2n)$ Work and Low Depth](https://dl.acm.org/doi/10.1145/3565557) in a shared memory machine. The referenced paper is a theory paper, and we aim for filling up the gap between theory and implementations. Specifically, we implement parallel randomized edge contractions using parallel tree contractions and parallel batched updates and queries. Although our implementation indeed scales in a shared memory machine, it's way slower than a seqential counterpart due to the high overhead stemming from the data structure used in parallel randomized edge contractions. Also, our implementation is outperformed by existing parallel deterministic solvers, which raises the question whether this parallel randomized algorithm is useful in practice. Please refer to [report.pdf](report.pdf) for more details of this project.

## Dependencies
- OpenMP
- GBBS (gbbs) as a submodule, use `cd gbbs && git apply ../gbbs.patch && cd ..` to remove debug information and fix issues in templates and atomic operations
- ParlayLib (parlaylib) as a submodule, use `cd parlaylib && git apply ../parlaylib.patch && cd ..` to fix issues in the scheduler
- PBBS (pbbsbench) as a submodule

## Compliation
```
$ make
```

## Usage
The binary is located at `./build/min-cut`.
```
$ ./build/min-cut -h

Usage: ./build/min-cut [options]

Options:
    -g <Source Graph File Path>                 Path of file with source graph in edge list format.
    -k <Number of rounds of edge contractions>  Default to 1, check the referenced paper to see how many rounds are needed to achieve the desired approximation w.h.p.
    -s <Seed>                                   Default to 42, the seed used for randomized algorithms.
    -p <Paradigm>                               Choose among parrec (default), seqrec-uf, and seqrec-rc. Check the report for more information.
    -h                                          Help.
```
One can use OpenMP environment variables and numactl as needed. An example run follows:
```
OMP_PLACES=cores numactl -i all ./build/min-cut -g examples/toy.txt -k 10 -p parrec
```
# Graph Format
We use the edge list format as follows, with an example given in [examples/toy.txt](examples/toy.txt):
```
# There can be comments at the top of the file, where each comment starts with '#'.
<source of edge 0> <destination of edge 0> <weight of edge 0>
<source of edge 1> <destination of edge 1> <weight of edge 1>
...
<source of edge m - 1> <destination of edge m - 1> <weight of edge m - 1>
```