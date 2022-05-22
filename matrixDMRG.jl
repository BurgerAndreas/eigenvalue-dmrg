#=

For debugging run these commands in the Julia REPL:
ITensors.enable_debug_checks()
ITensors.disable_debug_checks()
=#

#import Pkg; Pkg.add("stdlib")

using ITensors
using SparseArrays


function compute_eigenvalues(matrix; N = 20, numValues = 1)
    # (1) Define the matrix 
    # Define the vector space (2^size x 2^size)
    sites = siteinds("S=1",N)

    # Matrix in terms of Pauli matrices
    os = matrixToOpSum(matrix)
    # Todo placeholder 
    os = OpSum()
    for j=1:N-1
        os += -1,"Sz",j,"Sz",j+1
    end
    for j=1:N
        os += -0.5,"Sx",j;
    end
    H = MPO(os,sites)

    # (2) Configure DMRG
    sweeps = Sweeps(30)
    maxdim!(sweeps,10,10,10,20,20,40,80,100,200,200)
    cutoff!(sweeps,1E-8)
    noise!(sweeps,1E-6)
    weight = 20

    # (3) Run DMRG: Smallest Eigenvalue
    vector0_init = randomMPS(sites,linkdims=2)
    value0,vector0 = dmrg(H,vector0_init,sweeps)
    values = [value0]
    vectors = [vector0]

    # (4) Next smallest Eigenvalues
    for x=1:numValues
        vectorx_init = randomMPS(sites,linkdims=2)
        valuex,vectorx = dmrg(H,vectors,vectorx_init,sweeps; weight)
        push!(values, valuex)
        push!(vectors, vectorx)
    end
   
    return values
end


# Compressed Sparse Column (CSC) Sparse Matrix
function randomSparseMatrix(; symmetric = false, size = 20, sparsity = 0.1)
    if symmetric
        return Symmetric(sprand(size,size,sparsity))
    else
        return sprand(size,size,sparsity)
    end
end


# Todo placeholder
function matrixToOpSum(matrix)
    os = OpSum()

    # Row indices, column indices, values
    entries = findnz(matrix)

    # Loop over matrix entries
    for entry in eachindex(SparseArrays.nnz(matrix))
        # Add Paulis to OpSum
        row = entries[1][entry]
        col = entries[2][entry]
        value = entries[3][entry]
        
        println(row, " ", col, " ", value)
    end

    return os
end


# Program
let
    # Matrix
    N = 3
    size = 2^N
    a = randomSparseMatrix(size = size)
    println("matrix type", typeof(a))
    println(a)

    os = matrixToOpSum(a)

    # DMRG
    #values = compute_eigenvalues(matrix; N = 20, numValues = 2)
    #println("values", values)
end