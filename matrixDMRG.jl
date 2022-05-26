#=

For debugging run these commands in the Julia REPL:
ITensors.enable_debug_checks()
ITensors.disable_debug_checks()
=#

#import Pkg; Pkg.add("stdlib")

using ITensors
using SparseArrays
using LinearAlgebra # eigen(), eigvals()


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
        row_binary = reverse(digits(entries[1][entry]-1, base=2, pad=N))
        col_binary = reverse(digits(entries[2][entry]-1, base=2, pad=N))
        value = entries[3][entry]

        # Replace each row-column-digit-pair with a 2x2 matrix 
        # (made out of Pauli matrices)
        # (Number,String,Int,String,Int,...)
        os_vector1 = []
        os_vector2 = []
        push!(os_vector1, value)
        push!(os_vector2, value)
        # loop over digit-pairs
        for pair in 1:N
            if row_binary[pair] == 0
                if col_binary[pair] == 0
                    # 00 = "0.5(Id+2*Z)=Pu"
                    push!(os_vector1, "Id", pair) # operator
                    os_vector1[1] *= 0.5 # factor
                    push!(os_vector2, "Sz", pair) # operator
                else
                    # 01 = "0.5(X+iY)=P+=S+"
                    push!(os_vector1, "S+", pair) # operator
                    push!(os_vector2, "S+", pair) # operator
                end
            else
                if col_binary[pair] == 0
                    # 10 = "0.5(X-iY)=P-=S-"
                    push!(os_vector1, "S-", pair) # operator
                    push!(os_vector2, "S-", pair) # operator
                else
                    # 11 = "0.5(Id-2*Z)=Pd"
                    push!(os_vector1, "Id", pair) # operator
                    os_vector1[1] *= 0.5 # factor
                    push!(os_vector2, "Sz", pair) # operator
                    os_vector1[1] *= -1 # factor
                end
            end
        end
        # Add Pauli matrices to OpSum
        # OpSum only accepts (immutable) tuples
        os += tuple(os_vector1...)
        os += tuple(os_vector2...)

    end # loop over matrix entries

    return os
end

function show_operators()
    # https://itensor.github.io/ITensors.jl/dev/SiteType.html
    s = siteind("S=1/2")
    println("\nOperators:")
    @show op("Id", s);
    #@show op("Sx", s);
    #@show op("Sy", s);
    @show op("Sz", s);
    @show op("S+", s);
    @show op("S-", s);
end

function compute_eigenvalues(matrix, N; eigenvalues = 1)
    println("\n-----------------------------------")
    # (1) Define the matrix 
    # Define the vector space (2^size x 2^size)
    sites = siteinds("S=1/2",N)

    # Matrix in terms of Pauli matrices
    os = matrixToOpSum(matrix, N)
    H = MPO(os,sites)

    # (2) Configure DMRG
    sweeps = Sweeps(30)
    maxdim!(sweeps,10,10,10,20,20,40,80,100,200,200)
    cutoff!(sweeps,1E-8)
    noise!(sweeps,1E-6)
    weight = 20
    outputlevel = 0

    # (3) Run DMRG: Smallest Eigenvalue
    vector0_init = randomMPS(sites,linkdims=2)
    value0,vector0 = dmrg(H,vector0_init,sweeps;outputlevel)
    values = [value0]
    vectors = [vector0]

    # (4) Next smallest Eigenvalues
    for x=1:(eigenvalues-1)
        vectorx_init = randomMPS(sites,linkdims=2)
        valuex,vectorx = dmrg(H,vectors,vectorx_init,sweeps;weight,outputlevel)
        push!(values, valuex)
        push!(vectors, vectorx)
    end

    # Output
    #@show os
    println("Number of Pauli Operators: ", length(os))
   
    return values
end


# Program
let
    #= Define Matrix =#
    # Vector space: (2^N)x(2^N)
    N = 5 
    size = 2^N
    # Sparse Matrix in CSC format
    a = randomSparseMatrix(size = size)

    #= Run DMRG =#
    # DMRG
    @time values = compute_eigenvalues(a, N; eigenvalues=3)


    #= Compare with Eigensolver =#
    # https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/
    # vals_ex, vec_ex = eigen(a)
    @time exact_values = eigvals(Matrix(a))
    exact_values = sort(real(exact_values))

    # Output
    #println(a)
    println("\nEigenvalues DMRG ", values)
    println("Eigenvalues exact ", exact_values[1:3])

    #= Test components =#
    #show_operators()
    #@time os = matrixToOpSum(a, N)
    #println(os)

end