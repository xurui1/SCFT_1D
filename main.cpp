
#include "global.h"
#include "parameters.h"
#include "omega.h"
#include "TDMA.h"
#include "solvediffeq.h"
#include "vol.h"
#include "phi.h"
#include "Q_partition.h"
#include "Conc.h"
#include "Incomp.h"
#include "output.h"
#include "fE.h"
#include "FreeEnergy.h"
#include "homogfE.h"



int main( ){
    
    double **w;
    double *eta;
    double **phi;
    double *chi;
    double *f;
    double *mu;
    double ds;
    int *Ns;
    double dr;
    double volume;
    double **chiMatrix;
    double fE_hom;
    
    //Allocate memory
    w=create_2d_double_array(ChainType,Nr,"w");          //Auxiliary potential fields
    eta=create_1d_double_array(Nr,"eta");                //Incompressibility field
    phi=create_2d_double_array(ChainType,Nr,"phi");      //Concentration fields
    chi=create_1d_double_array(ChainType,"chi");            //Interaction parameters
    f=create_1d_double_array(ChainType,"f");                //Chain fractions
    Ns=create_1d_integer_array(ChainType, "Ns");            //Chain lengths
    mu=create_1d_double_array(2, "mu");                     //Chemical potentials
    chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");
    
    //Initial time for random number generator
    long iseed;
    time_t t;
    iseed=time(&t);
    srand48(iseed);
    
    //Set parameters
    parameters(chi,f,&ds,Ns,&dr,chiMatrix,mu);
    volume=vol(dr);
    
    //Calculate homogeneous free energy
    fE_hom=homogfE(mu,chiMatrix,f);
    
    //Set up initial omega field
    omega(w);
    
    //SCFT
    FreeEnergy(w,phi,eta,Ns,ds,chi,dr,chiMatrix,mu,fE_hom,volume);
    
    //Destroy memory allocations------------
    destroy_2d_double_array(w);
    destroy_1d_double_array(eta);
    destroy_2d_double_array(phi);
    destroy_1d_double_array(chi);
    destroy_1d_integer_array(Ns);
    destroy_1d_double_array(f);
    destroy_2d_double_array(chiMatrix);
    //-------------------------------------
    
    return 0;
}
