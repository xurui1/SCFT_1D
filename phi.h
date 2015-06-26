void phi_total(double *phiA, double *phiB, double *phiC, double dr, double volume){
    
    //Here I am calculating the total concentration of each species using a trapezoidal (?) rule. I removed the trapezoidal component because I don't think it should matter too much if we are using the Neumann boundary condition. Alternatively, I think I would have to modify the incompressibility condition.
    
    double phiA_tot,phiB_tot,phiC_tot;
    phiA_tot=0.0;
    phiB_tot=0.0;
    phiC_tot=0.0;
    double phi_tot=0.0;
    int i;
    
    
    //ends
    phiA_tot+=0.5*phiA[0]*dV(0,dr);
    phiB_tot+=0.5*phiB[0]*dV(0,dr);
    phiC_tot+=0.5*phiC[0]*dV(0,dr);
    
    phiA_tot+=0.5*phiA[(int)Nr-1]*dV(Nr-1,dr);
    phiB_tot+=0.5*phiB[(int)Nr-1]*dV(Nr-1,dr);
    phiC_tot+=0.5*phiC[(int)Nr-1]*dV(Nr-1,dr);

 
    for (i=1;i<(int)Nr-1;i++){
            phiA_tot+=phiA[i]*dV(i,dr);
            phiB_tot+=phiB[i]*dV(i,dr);;
            phiC_tot+=phiC[i]*dV(i,dr);;
    }
    
    //normalize
    phiA_tot/=volume;
    phiB_tot/=volume;
    phiC_tot/=volume;
    phi_tot=phiA_tot+phiB_tot+phiC_tot;
    
    cout<<"phiA: "<<phiA_tot<<" phiB: "<<phiB_tot<<" phiC: "<<phiC_tot<<" total: "<<phi_tot<<endl;
    
    
};


void phi_calc(double **phi,double **qA,double **qdagA,double **qB,double **qdagB,double **qC,int *Ns,double*mu,double ds){
    int i,s;
    
    for(i=0;i<Nr;i++){
            //Empty array elements
            phi[0][i]=0.0;
            phi[1][i]=0.0;
            phi[2][i]=0.0;
            
            //phiA integration
            for(s=0;s<(int)Ns[0]+1;s++){
                if(s==0 || s==(int)Ns[0]){
                    phi[0][i]+=0.5*qA[i][s]*qdagA[i][Ns[0]-s]*ds;
                }
                else{
                    phi[0][i]+=qA[i][s]*qdagA[i][Ns[0]-s]*ds;
                }
            }
            
            //phiB integration
            for(s=0;s<(int)Ns[1]+1;s++){
                if(s==0 || s==(int)Ns[1]){
                    phi[1][i]+=0.5*qB[i][s]*qdagB[i][Ns[1]-s]*ds;
                }
                else{
                    phi[1][i]+=qB[i][s]*qdagB[i][Ns[1]-s]*ds;
                }
            }
            
            //phiC integration
            for(s=0;s<(int)Ns[2]+1;s++){
                if(s==0 || s==(int)Ns[2]){
                    phi[2][i]+=0.5*qC[i][s]*qC[i][Ns[2]-s]*ds;
                }
                else{
                    phi[2][i]+=qC[i][s]*qC[i][Ns[2]-s]*ds;
                }
            }
            
            //Grand canonical relation
            phi[0][i]=exp(mu[0])*phi[0][i];
            phi[1][i]=exp(mu[0])*phi[1][i];
            phi[2][i]=exp((mu[1])*kappa)*phi[2][i]*(1.0/kappa);
    }
    
    
};
