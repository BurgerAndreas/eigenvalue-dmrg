#include "itensor/all.h"

using namespace itensor;

/*
 Andreas Burger, Heitor Casagrande, May 2022
 
 Find the smallest Eigenvalues of a hermitian matrix.
 Using DMRG and the ITensor library.
 (1) Define the matrix
 (2) Configure DMRG
 (3) Run DMRG: Smallest Eigenvalue
 
 (4) 2nd smallest Eigenvalue
 (5) 3rd smallest Eigenvalue
 (6) 4th smallest Eigenvalue
 (7) ...
 */


int main()
    {
    /* (1) Define the matrix */
    // Matrix size
    int size = 32;

    // Define the vector space (2^size x 2^size)
    auto sites = SpinOne(size,{"ConserveQNs=",false});

    // Matrix parameters
    Real h = 1.0;
    Real J = 1.0;

    // Matrix in terms of Pauli matrices
    auto ampo = AutoMPO(sites);
    for(int j = 1; j < size; ++j){
        ampo += -J,"Sz",j,"Sz",j+1;
    }
    for(int j = 1; j <= size; ++j){
        ampo += -h,"Sx",j;
    }
    auto H = toMPO(ampo);

    
    /* (2) Configure DMRG */
    auto sweeps = Sweeps(30);
    sweeps.maxdim() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    println(sweeps);

    
    /* (3) Run DMRG: Smallest Eigenvalue */
    // Run DMRG, return Eigenvalue and Eigenvector
    auto [value0, vector0] = dmrg(H,randomMPS(size),sweeps,{"Quiet=",true});

    // Print the smallest Eigenvalue
    println("\n----------------------\n");
    printfln("Smallest Eigenvalue = %.10f",value0);


    /* (4) 2nd smallest Eigenvalue */
    // transform Eigenvector0 into an MPS that the DMRG can read
    auto vector0_mps = std::vector<MPS>(1);
    vector0_mps.at(0) = vector0;

    // Run DMRG
    auto [value1, vector1] = dmrg(H,vector0_mps,randomMPS(sites),sweeps,{"Quiet=",true,"Weight=",20.0});

    // Print the 2nd smallest Eigenvalue
    printfln("2nd smallest Eigenvalue = %.10f",value1);
    // Test: the overlap between the Eigenvectors should be very close to zero
    printfln(" Overlap of Eigenvectors: vector0*vector1 = %.2E",inner(vector0,vector1));
    
    
    /* (5) 3rd smallest Eigenvalue */
    // transform Eigenvector0 and Eigenvector1 into an MPS that the DMRG can read
    auto vector01_mps = std::vector<MPS>(2);
    vector01_mps.at(0) = vector0;
    vector01_mps.at(1) = vector1;

    // Run DMRG
    auto [value2, vector2] = dmrg(H,vector01_mps,randomMPS(sites),sweeps,{"Quiet=",true,"Weight=",20.0});

    // Print the 3rd smallest Eigenvalue
    printfln("3rd smallest Eigenvalue = %.10f",value2);
    // Test: the overlap between the Eigenvectors should be very close to zero
    printfln(" Overlap of Eigenvectors: vector1*vector2 = %.2E",inner(vector1,vector2));
    printfln(" Overlap of Eigenvectors: vector0*vector2 = %.2E",inner(vector0,vector2));
    
    
    /* (6) 4th smallest Eigenvalue */
    // ...
    

    return 0;
    }
