# Efficient Maximum s-bundle Search Via Local Vertex Connectivity

This repository implements the maximum *s*-bundle computation algorithm **SymBD** proposed in our SIGMOD 2025 paper. If you are using the code, please cite our paper.


## Compile the code

```sh
$ make clean
$ make
```
It generates an executable program "SymBD", which corresponds to our SymBD algorithm (**Algorithm 3** in the paper).

## Get datasets
The three dataset collections used in our paper can be downloaded from the links below.

1. 2nd DIMACS graphs collection: http://archive.dimacs.rutgers.edu/pub/challenge/graph/benchmarks/clique/

2. 10th DIMACS graphs collection: https://networkrepository.com/dimacs10.php

3. Real-world graphs collection: https://networkrepository.com/index.php


## Data format
Our program only accepts a **binary input data file** (i.e., with a **.bin** file extension). If the data is not in binary form, you need to use our provided data transformation tool `translate.cpp` which is contained in the directory `data_transform_tool` to transform the data. Additionally, we have included an **original data file** `bio-celegans.mtx` within the same file directory, which can assist in understanding the native format of the graph data to be read.

We provide an example of converting a data file with a '.mtx' suffix into the binary data below. 

    g++ -O3 -std=c++11 translate.cpp -o a.out -w
    ./a.out bio-celegans.mtx bio-celegans.bin

We remark that processing approaches for other types of graph data formats are similar. **If any issues related to data conversion arise during this process, please refer to the very simple and clear source code of our conversion function and make the necessary modifications.** Furthermore, to provide a clearer understanding of the data format prior to conversion to binary representation, we present a simple illustrative example as shown below (note that the two integers in the first line means **the number of vertices and edges**, respectively).

    4 5
    1 2
    1 4
    2 3
    2 4
    3 4
    
## Run the code
The usage procedure for the generated executable program is as follows.
```sh
$ ./SymBD -g {path_to_graph} -s {s_value}
```

To illustrate, an example of computing the exact maximum 6-bundle for our provided example graph 'bio-pdb1HYS.bin' is as follows
```sh
$ ./SymBD -g datasets/bio-pdb1HYS.bin -s 6
```
## Result Analysis
If you are using our provided example above, then you will see the following outputs (may be there are some slight differences):
```sh
**** SymBD (Release) build at 16:11:27 Jul 11 2024 ***
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
From the obtained empirical results, we can extract the following key information:

- Graph data characteristics:
  - File Name: `bio-pdb1HYS`
  - Vertices (n): 36,417
  - Edges (m): 2,154,174

- s-Bundle sizes bbtained via our three-stage heuristic approach (The details can be found in **Algorithm 4**, Section 5.2 of our paper.):
  - Stage 1: Degeneracy s-Bundle Size: 58
  - Stage 2: 1-hop Subgraph Degen s-Bundle Size: 60
  - Stage 3: 2-hop-Degen s-Bundle Size: 63

- Performance metrics:
  - Preprocessing time: 1,286,236 microseconds
  - Generated branches: 518
  - Maximum 6-bundle size: 63
  - Total runtime: 1,362,248 microseconds

---

##  Parallel Version
We further provide a parallel version of the code for solving the maximum s-bundle problem, which is located in the ***parallel_version*** folder.

 To compile the project, you should navigate to the folder and then use the following commands:

```sh
$ g++ -std=c++20 -O3 -fopenmp -o SymBD_parallel main.cpp Graph.cpp sbundle_tool.cpp
```

If you wish to run the project that you have just compiled, you can use the graph data '*ia-wiki-user-edits-page.bin*' that we provide in the folder. The usage procedure for the generated executable program is as follows.
```sh
$ ./SymBD_parallel {path_to_graph} {s_value} {num_of_threads}
```

To illustrate, the following is an example of parallel computation of the exact maximum 6-bundle for our provided example graph '*ia-wiki-user-edits-page.bin*' using 8 threads:

```sh
$ ./SymBD_parallel ia-wiki-user-edits-page.bin 6 8
```

Note: **if you encounter any issues throughout the process, please do not hesitate to contact us. My contact information is: ly17369279121@163.com.**
