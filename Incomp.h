void Incomp(double *eta, double **phi, double *delphi){
    
    int     i;
    int     chain;
    double  ptot;
    
    ptot=0.0;
    
    for(i=0;i<Nr;i++){
        
        ptot=0.0;
        delphi[i]=0.0;
                
        for(chain=0;chain<ChainType;chain++){
            ptot+=phi[chain][i];
        }
                            
        delphi[i]=1.0-ptot;
        eta[i]-=delphi[i];
        
        if (fabs(delphi[i])>1e2){
            cout<<i<<" incomp: "<<delphi[i]<<endl;
            exit(EXIT_FAILURE);
        }
    }
}
