void FreeEnergy(double **w, double **phi, double *eta, int *Ns, double ds, double *chi, double dr, double **chiMatrix, double *mu, double fE_hom, double volume){
    
    
    double  currentfE, oldfE, deltafE;
    int     maxIter=10000;
    double precision=1e-6;          //convergence condition
    int     i,iter,chain,ii,jj;
    double  Q;
    double  fE_int, fES;            //interaction free energy and chain partition function fE
    double  epsilon, gamma;
    double  *delphi;
    double  **delW;
    double  **newW;
    double  deltaW;
    
    //Arrays for updating the omega fields
    delW=create_2d_double_array(ChainType,Nr,"delW");
    delphi=create_1d_double_array(Nr,"delphi");
    newW=create_2d_double_array(ChainType,Nr,"newW");
    
    currentfE=0.0;
    deltafE=0.0;
    
    epsilon=0.05;
    gamma=0.05;
    
    iter=0;
    std::ofstream outputFile("./results/fE.dat");
    for (iter=0;iter<maxIter;iter++){
        
        fE_int=0.0;
        fES=0.0;
        deltaW=0.0;

        
        Q=Conc(phi,w,Ns,ds,dr,mu,volume);          //Calculate Chain partition function for both AB and C
        
        
        Incomp(eta,phi,delphi);              //Enforce incompressibility condition
        output(dr,phi,w);                   //Output some data to file
    
        
        //clear omega field update
        for(i=0;i<Nr;i++){
            for(ii=0;ii<ChainType;ii++){
                newW[ii][i]=0.0;
            }
        }
        
        //Calculate components for new field and interaction free energies
        for(i=0;i<Nr;i++){
            for(ii=0;ii<ChainType;ii++){
                for(jj=0;jj<ChainType;jj++){
                    newW[ii][i]+=((chiMatrix[ii][jj]*phi[jj][i])+eta[i]);
                }
                delW[ii][i]=newW[ii][i]-w[ii][i];
                deltaW+=fabs(delW[ii][i])*dV(i,dr);
                }
        }
        fE_int=fE(newW,phi,chiMatrix,dr,volume);
        
        //Normalize by box size
        deltaW/=volume;
        
       
        //Update free energy
        fES=Q;
        oldfE=currentfE;
        currentfE=-fES+fE_int;
        
        if (fabs(currentfE)>1e3){
            cout<<"fE too large"<<endl;
            exit(EXIT_FAILURE);
        }
        deltafE=fabs(currentfE-oldfE);
        
        //Print free energy, difference in free energy, change in omega field to screen
        std::cout<<iter<<" fE:"<<currentfE<< " dfE:"<<currentfE-fE_hom<<" " << deltaW<<" "<<fE_hom<<std::endl;
        outputFile << iter << " " << currentfE<< " " << currentfE-fE_hom<<" "<< deltaW<<std::endl;
        
        //Update omega field. Can we add Anderson mixing for this system?
        for(i=0;i<Nr;i++){
            for(chain=0;chain<ChainType;chain++){
                w[chain][i]+=(gamma*delW[chain][i]-epsilon*delphi[i]);
            }
        }
    
        if (deltafE<precision && deltaW<2.0){break;} //Convergence condition
        
    }
    
    outputFile.close();
    
    destroy_1d_double_array(delphi);
    destroy_2d_double_array(delW);
    destroy_2d_double_array(newW);
    
}
