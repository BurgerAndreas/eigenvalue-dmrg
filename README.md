# eigenvalue-dmrg
Eigenvalues of large matrices using DMRG and ITensor.

ITensor - Software Library for Tensor Network Calculations <br />
by Matthew Fishman and Steven R. White and E. Miles Stoudenmire <br />
arxiv.org/abs/2007.14822 <br />
http://itensor.org/

In collaboration with: <br />
Heitor Peres Casagrande https://github.com/heitorc7

### DMRG for any matrix
matrixDMRG.jl

### Examples
simpleDMRG.jl - Single smallest Eigenvalue of a Hamiltonian using DMRG <br />
eigenvaluesDMRG.jl - First few smallest Eigenvalue of a Hamiltonian using DMRG <br />
eigenvaluesDMRG.cc - Same, but using ITensor C++ version <br />

# Installation

### Julia
https://julialang.org/downloads/

### ITensor Julia version
https://itensor.github.io/ITensors.jl/stable/getting_started/Installing.html <br />
Open Julia (REPL) <br />
julia> ] <br />
pkg> add ITensors <br />
press backspace <br />
julia> using ITensors; ITensors.compile() <br />

### Optional: ITensor C++ version
(0) Download and build ITensor <br />
http://itensor.org/docs.cgi?vers=cppv3&page=install <br />
$ git clone https://github.com/ITensor/ITensor itensor <br />
$ cd itensor <br />
$ cp options.mk.sample options.mk <br />
Edit options.mk <br />
$ make 

(1) Change path to ITensor <br />
Makefile -> change line 3

(2) Build and execute <br />
$ cd eigenvalue-dmrg/Examples <br />
$ make <br />
$ ./eigenvaluesDMRG
