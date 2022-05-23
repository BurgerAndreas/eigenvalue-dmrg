#=

For debugging run these commands in the Julia REPL:
ITensors.enable_debug_checks()
ITensors.disable_debug_checks()
=#

#import Pkg; Pkg.add("stdlib")

using ITensors
using SparseArrays


function compute_eigenvalues(matrix; N = 20, eigenvalues = 1)
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
    for x=1:eigenvalues
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


# Todo: complex
# Input: (2^N)x(2^N) matrix
# Output: 
# Each entry is replaced by a tensorproduct of N 2x2 matrices
# Each 2x2 matrix is a linear combination of Pauli matrices and the Idendity
function matrixToOpSum(matrix, N)
    os = OpSum()
    # Todo: proper OpSum not implemented
    os = []

    # Row indices, column indices, values
    entries = findnz(matrix)

    # Loop over matrix entries
    for entry in 1:SparseArrays.nnz(matrix)
        # Write row and column index in binary digits (length N)
        row_binary = digits(entries[1][entry], base=2, pad=N)
        col_binary = digits(entries[2][entry], base=2, pad=N)
        value = entries[3][entry]

        # Replace each row-column-digit-pair with a 2x2 matrix 
        # (made out of Pauli matrices)
        # Add Pauli matrices to OpSum
        # Todo: probably not easily possible
        # Decompose tensorproduct of N matrices as a sum of tensorproducts of 2 matrices?
        for pair in 1:N
            if row_binary[pair] == 0
                if col_binary[pair] == 0
                    # 00
                    push!(os, "0.5(Id+Z)")
                else
                    # 01
                    push!(os, "0.5(X+iY)=P+")
                end
            else
                if col_binary[pair] == 0
                    # 10
                    push!(os, "0.5(X+iY)=P-")
                else
                    # 11
                    push!(os, "0.5(Id-Z)")
                end
            end
        end

    end

    return os
end


# Program
let
    # Matrix
    # Vector space: (2^N)x(2^N)
    N = 3 
    size = 2^N
    # Sparse Matrix in CSC format
    a = randomSparseMatrix(size = size)
    println("matrix type", typeof(a))
    println(a)

    os = matrixToOpSum(a, N)
    println(os)

    # DMRG
    #values = compute_eigenvalues(a; N = 20, eigenvalues = 2)
    #println("values", values)
end