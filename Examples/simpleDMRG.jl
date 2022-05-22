# Example from
# https://itensor.github.io/ITensors.jl/stable/

using ITensors
let
  # Create 100 spin-one indices
  N = 32
  sites = siteinds("S=1",N)

  # Input operator terms which define
  # a Hamiltonian matrix, and convert
  # these terms to an MPO tensor network
  # (here we make the 1D Heisenberg model)
  ampo = OpSum()
  for j=1:N-1
    ampo += "Sz",j,"Sz",j+1
    ampo += 0.5,"S+",j,"S-",j+1
    ampo += 0.5,"S-",j,"S+",j+1
  end
  H = MPO(ampo,sites)

  # Create an initial random matrix product state
  psi0 = randomMPS(sites)

  # Plan to do 5 passes or 'sweeps' of DMRG,
  # setting maximum MPS internal dimensions
  # for each sweep and maximum truncation cutoff
  # used when adapting internal dimensions:
  sweeps = Sweeps(5)
  setmaxdim!(sweeps, 10,20,100,100,200)
  setcutoff!(sweeps, 1E-10)
  @show sweeps

  # Run the DMRG algorithm, returning energy
  # (dominant eigenvalue) and optimized MPS
  energy, psi = dmrg(H,psi0, sweeps)
  println("Final energy = $energy")

  nothing
end