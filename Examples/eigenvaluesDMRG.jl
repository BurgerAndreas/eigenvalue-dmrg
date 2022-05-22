# Example from
# https://itensor.github.io/ITensors.jl/stable/examples/DMRG.html

using ITensors

function compute_eigenvalues(; N = 20, h = 1.0, J = 0.5)
    # (1) Define the matrix 
    #Define the vector space (2^size x 2^size)
    sites = siteinds("S=1",N)

    # Matrix in terms of Pauli matrices
    # Use the OpSum feature to create the
    # transverse field Ising model
    os = OpSum()
    for j=1:N-1
        os += -J,"Sz",j,"Sz",j+1
    end
    for j=1:N
        os += -h,"Sx",j;
    end
    H = MPO(os,sites)

    # (2) Configure DMRG
    # Make sure to do lots of sweeps
    # when finding excited states
    sweeps = Sweeps(30)
    maxdim!(sweeps,10,10,10,20,20,40,80,100,200,200)
    cutoff!(sweeps,1E-8)
    noise!(sweeps,1E-6)
    weight = 20*h

    # (3) Run DMRG: Smallest Eigenvalue
    # Run DMRG, return Eigenvalue and Eigenvector
    vector0_init = randomMPS(sites,linkdims=2)
    value0,vector0 = dmrg(H,vector0_init,sweeps)
    println("Smallest Eigenvalue =", value0)

    # (4) 2nd smallest Eigenvalue
    vector1_init = randomMPS(sites,linkdims=2)
    value1,vector1 = dmrg(H,[vector0],vector1_init,sweeps; weight)
    println("2nd smallest Eigenvalue =", value1)

    # Test: the overlap between the Eigenvectors should be very close to zero
    @show inner(vector1,vector0)
    println()

    # (5) 3rd smallest Eigenvalue
    vector2_init = randomMPS(sites,linkdims=2)
    value2,vector2 = dmrg(H,[vector0,vector1],vector2_init,sweeps;weight)
    println("2nd smallest Eigenvalue =", value2)

    # Test: the overlap between the Eigenvectors should be very close to zero
    @show inner(vector2,vector0)
    @show inner(vector2,vector1)

    return
end

compute_eigenvalues()