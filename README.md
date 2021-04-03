# The Hierarchical Poincare-Steklov (HPS) Method for Elliptic Partial Differential Equations

A work in progress!

Submitted April 2nd, 2021 as part of my Comprehensive Exam for the PhD in Computing - Computational Sciences and Engineering at Boise State University

## Created by: Damyn Chipman (aka [camperD](https://github.com/camperD/))

---

- BSU School of Computing
- Computational Sciences and Engineering
- DamynChipman@u.boisestate.edu

---

### Introduction

This repo contains progress and current implementation of the HPS method for a fast and direct elliptic PDE solver to be used on adaptive meshes.

This is ongoing and part of my dissertation research. Updates will come as implemented.

#### Rough Timeline

* March 1-5 2021 : SIAM CSE presentation on progress of method and results.
* Spring/Summer 2021 : Parallel implementation of HPS method, incorporation with `p4est` and `ForestClaw`.
* More to come!

### Installing

Install using CMake. Clone or download this repo:

```bash
$ git clone https://github.com/camperD/HPS
```

Then perform the following standard CMake commands:

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make test
```

### Demo Usage

Running `make` will generate three executables in `${build_directory}/bin`: `convergence`, `4to1_merge`, and `HPS`. Each provides further testing and demonstrates the HPS method.

#### Convergence Analysis

Test for patch convergence with:

```bash
$ convergence <Poisson Problem ID>: 0, 1, 2, 3, 4>
```

The program will use the patch solver on increasing number of cells and check the error with a known solution. The current patch solver has 2nd order convergence (implements FISHPACK90).

#### 4-to-1 Merge Testing

The `4to1_merge` program will create 4 leaf patches and merge horizontally and vertically to form a parent patch. It will then check if the operators created in the merge process are equal to the operators created on a refined parent level grid.

```bash
$ 4to1_merge <Poisson Problem ID: 0, 1, 2, 3, 4>
```

#### HPS Build and Solve

This is the primary executable. `HPS` will create a 3-level tree (4 leaf patches) and then perform the build and solve stages of the HPS method. It will compute the error compared to the patch solver implemented.

```bash
$ HPS <Number of leaf cells per side> <Poisson Problem ID: 0, 1, 2, 3, 4> <output files prefix>
```

Note: Currently implemented for Laplace's equation (IDs: 0, 1, and 2). Poisson equation implementation is being worked on currently (April 2nd, 2021).

### References and Acknowledgements

Dr. Donna Calhoun - Research Advisor and Bug Finder Extraordinaire

Martinsson P.G., 2015, 'The Hierarchical Poincare-Steklov (HPS) solver for elliptic PDEs: A tutorial'.

Martinsson, Per-Gunnar. *Fast Direct Solvers for Elliptic PDEs*. Philidelphia, PA, SIAM, 2020.
