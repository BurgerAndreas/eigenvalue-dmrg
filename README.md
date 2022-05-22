# eigenvalue-dmrg
Eigenvalues of large matrices using DMRG and ITensor.

ITensor - Software Library for Tensor Network Calculations
by Matthew Fishman and Steven R. White and E. Miles Stoudenmire
arxiv.org/abs/2007.14822
http://itensor.org/

In collaboration with:
Heitor Peres Casagrande https://github.com/heitorc7

### DMRG for any matrix
matrixDMRG.jl

### Examples
simpleDMRG.jl - Single smallest Eigenvalue of a Hamiltonian using DMRG
eigenvaluesDMRG.jl - First few smallest Eigenvalue of a Hamiltonian using DMRG
eigenvaluesDMRG.cc - Same, but using ITensor C++ version

# Installation

### Julia
https://julialang.org/downloads/

### ITensor Julia version
https://itensor.github.io/ITensors.jl/stable/getting_started/Installing.html
Open Julia (REPL)
julia> ]
pkg> add ITensors
press backspace
julia> using ITensors; ITensors.compile()

### Optional: ITensor C++ version
(0) Download and build ITensor
http://itensor.org/docs.cgi?vers=cppv3&page=install
$ git clone https://github.com/ITensor/ITensor itensor
$ cd itensor
$ cp options.mk.sample options.mk
Edit options.mk
$ make

(1) Change path to ITensor
Makefile -> change line 3

(2) Build and execute
$ cd eigenvalue-dmrg/Examples
$ make
$ ./eigenvaluesDMRG
