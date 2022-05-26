# eigenvalue-dmrg
Eigenvalues of large matrices using DMRG and ITensor.

ITensor - Software Library for Tensor Network Calculations <br />
by Matthew Fishman and Steven R. White and E. Miles Stoudenmire <br />
arxiv.org/abs/2007.14822 <br />
http://itensor.org/


## DMRG for any matrix
matrixDMRG.jl

### Simpler examples
simpleDMRG.jl - Single smallest Eigenvalue of a Hamiltonian using DMRG <br />
eigenvaluesDMRG.jl - First few smallest Eigenvalue of a Hamiltonian using DMRG <br />
eigenvaluesDMRG.cc - First few smallest Eigenvalue of a Hamiltonian using DMRG, but using ITensor C++ version <br />

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
http://itensor.org/docs.cgi?vers=cppv3&page=install <br />
$ git clone https://github.com/ITensor/ITensor itensor <br />
$ cd itensor <br />
$ cp options.mk.sample options.mk <br />
Edit options.mk <br />
$ make <br />
Makefile -> Change path to ITensor in line 3 <br />
$ cd eigenvalue-dmrg/Examples <br />
$ make <br />
$ ./eigenvaluesDMRG

# Explanation
We rewrite a $2^N \times 2^N$ input matrix as a sum of tensorproducts of N matrices of size $2 \times 2$ <br />
$$ M_{2^N \times 2^N} \ =\ \sum_{entries}\ \otimes_N\ m_{2 \times 2}$$  
Where the $m_{2 \times 2}$ are build out of $Id, \sigma^x, \sigma^y, \sigma^z$ <br />

(1) Take each non-zero entry <br />
(2) Rewrite row and column numbers into binary digits <br />
(3) Each row-column-digit-pair corresponds to one $2 \times 2$ matrix <br />
(4) Tensorprodcut each matrix of a row-column-digit-pair ($N$ matrices per entry) <br />
Explained e.g. here: https://arxiv.org/abs/1909.12847

Problem:<br /> 
Probably not scalable to large N
