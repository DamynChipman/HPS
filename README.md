# The Hierarchical Poincare-Steklov (HPS) Method for Elliptic Partial Differential Equations

A work in progress!

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

* December 2020 : Initial public release of build and solve stages of HPS method for Poisson's equation.
* January 2021 - March 2021 : Integration with `p4est` package and adaptive mesh implementation.
* March 1-5 2021 : SIAM CSE presentation on progress of method and results.
* March 2021 - May 2021 : Parallel implementation of HPS method.
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

### References and Acknowledgements

Dr. Donna Calhoun - Research Advisor and Bug Finder Extraordinaire

Martinsson P.G., 2015, 'The Hierarchical Poincare-Steklov (HPS) solver for elliptic PDEs: A tutorial'.

Martinsson, Per-Gunnar. *Fast Direct Solvers for Elliptic PDEs*. Philidelphia, PA, SIAM, 2020.
