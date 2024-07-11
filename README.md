# Efficient Maximum s-bundle search via local vertex connectivity

This repository implements the maximum s-bundle computation algorithm **SymBK_H** proposed in our paper. If you are using the code, please cite our paper.

<!-- This repository is the same as https:XXX -->

## Compile the code

```sh
$ make clean
$ make
```
It generates an executable "SymBK_H", which corresponds to our SymBK_H algorithm.

## Get datasets
The data used in our paper are mainly from the links below.

1. 2nd DIMACS graphs collection: http://archive.dimacs.rutgers.edu/pub/challenge/graph/benchmarks/clique/

2. 10th DIMACS graphs collection: https://networkrepository.com/dimacs10.php

3. Real-world graphs collection: https://networkrepository.com/index.php


## Data format
Note that this program only accepts a **binary input data file** (with a .bin file extension). If the file does not end with `.bin`, you need to use the data format conversion function provided converter `translate.cpp` which is contained in the directory `data_transform_tool`.

Example for transform .mtx file to .bin file.

    g++ -g -o translate translate.cpp
    ./translate bio-celegans.mtx bio-celegans.bin

## Run the code
```sh
$ ./SymBK_H -g {path_to_graph} -s {s_value}
```

An example of computing the exact maximum 6-bundle for the graph bio-pdb1HYS.bin is as follows
```sh
$ ./SymBK_H -g datasets/bio-pdb1HYS.bin -s 6
```
## Result Analysis
If you are using our provided example above, then you can see the following outputs (may be there are some slight differences):
```sh
**** SymBK_H (Release) build at 16:11:27 Jul 11 2024 ***
File: bio-pdb1HYS n= 36417 m= 2154174 s = 6
Graph init ok
Now we find a larger solution, and its size is 59
*** Degeneracy s-bundle size: 58, max_core: 74, UB: 80, Time: 32,257 (microseconds)
*** After core shrink: n = 33,929, m = 2,048,332 (undirected)
*** After core_truss_copruning: n = 3279, m = 127259 (undirected)
*** 1-hop subgraph degen s-bundle size: 60, Time: 3,119 (microseconds)
*** After core_truss_copruning: n = 3279, m = 127259 (undirected)
2hop-Degen find a larger solution of size 61
2hop-Degen find a larger solution of size 62
2hop-Degen find a larger solution of size 63
*** 2-hop-Degen sbundle size: 63, Time: 13,944 (microseconds)
*** After core_truss_copruning: n = 345, m = 11289 (undirected)
*** The time for preprocessing is 1,286,236 (microseconds)
search_cnt: 10, ave_density: 0.97234, min_density: 0.93862
Branch_nodes: 518
*** Search time: 75,990
	Maximum 6-bundle Size: 63, Total Time: 1,362,248 (microseconds)
```
From the aforementioned empirical results, **we can extract the following key information**: the detailed characteristics of the processed graph data file (File: bio-pdb1HYS n= 36417 m= 2154174), the size of the s-bundle obtained via the heuristic approach (Degeneracy s-bundle size: 58, 1-hop subgraph degen s-bundle size: 60, 2-hop-Degen sbundle size: 63), the time consumed by the preprocessing stage (The time for preprocessing is 1,286,236 (microseconds)), the number of branches generated (Branch_nodes: 518), the size of the maximum s-bundle (Maximum 6-bundle Size: 63), and the total runtime of the algorithm (Total Time: 1,362,248 (microseconds)).
