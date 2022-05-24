#=

For debugging run these commands in the Julia REPL:
ITensors.enable_debug_checks()
ITensors.disable_debug_checks()
=#

#import Pkg; Pkg.add("stdlib")

using ITensors
using SparseArrays


function compute_eigenvalues(matrix, N; eigenvalues = 1)
    # (1) Define the matrix 
    # Define the vector space (2^size x 2^size)
    sites = siteinds("S=1",N)

    # Matrix in terms of Pauli matrices
    os = matrixToOpSum(matrix, N)
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


# Input: (2^N)x(2^N) matrix
# Output: Sum of tensorproduct of N 2x2 matrices
# Each entry is replaced by a tensorproduct of N 2x2 matrices
# Each 2x2 matrix is a linear combination of Pauli matrices and the Idendity
# Todo: complex entries
function matrixToOpSum(matrix, N)
    # Output
    os = OpSum()

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
        # (Number,String,Int,String,Int,...)
        os_vector = [] # vector for S+, S-
        os_vector_id = [] # vector for Id
        os_vector_plus = [] # vector for +Z
        os_vector_minus = [] # vector for  -Z
        for pair in 1:N
            # values Todo: fix prefactors
            # https://github.com/ITensor/ITensors.jl/blob/72b0bd9b2c9c2f37d311a4e3958f20d08ad7c069/src/physics/site_types/spinone.jl
            if pair == 1
                push!(os_vector, value) # value
                push!(os_vector_id, value*0.5) # value
                push!(os_vector_plus, value*0.5) # value
                push!(os_vector_minus, -value*0.5) # value
            end
            if row_binary[pair] == 0
                if col_binary[pair] == 0
                    # 00 = "0.5(Id+Z)"
                    push!(os_vector_id, "Id") # operator
                    push!(os_vector_id, pair) # position
                    push!(os_vector_plus, "Sz") # operator
                    push!(os_vector_plus, pair) # position
                else
                    # 01 = "0.5(X+iY)=P+"
                    push!(os_vector, "S+") # operator
                    push!(os_vector, pair) # position
                end
            else
                if col_binary[pair] == 0
                    # 10 = "0.5(X+iY)=P-"
                    push!(os_vector, "S-") # operator
                    push!(os_vector, pair) # position
                else
                    # 11 = "0.5(Id-Z)"
                    push!(os_vector_id, "Id") # operator
                    push!(os_vector_id, pair) # position
                    push!(os_vector_minus, "Sz") # operator
                    push!(os_vector_minus, pair) # position
                end
            end
        end
        # Add Pauli matrices to OpSum
        # OpSum only accepts (immutable) tuples
        vectors = [os_vector, os_vector_id, os_vector_plus, os_vector_minus]
        for vec in 1:length(vectors)
            if length(vectors[vec]) != 1
                #println("tuple ", tuple(vectors[vec]...))
                os += tuple(vectors[vec]...)
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

    #@time os = matrixToOpSum(a, N)
    #println(os)

    # DMRG
    @time values = compute_eigenvalues(a, N; eigenvalues = 1)
    #println("values", values)
end