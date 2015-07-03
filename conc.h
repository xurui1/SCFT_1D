double Conc(double **phi,double **w,int *Ns,double ds,double dr, double *mu,double volume){
    
    int         i,s;
    double      Q;
    double      **qA,**qB,**qC,**qdagA,**qdagB;
    double      *qintA,*qintB,*qintC;
    
    
    
    //Forwards propagators
    qA=create_2d_double_array(Nr,Ns[0]+1,"qA");
    qB=create_2d_double_array(Nr,Ns[1]+1,"qB");
    qC=create_2d_double_array(Nr,Ns[2]+1,"qC");
    
    //Complementary propagators
    qdagA=create_2d_double_array(Nr,Ns[0]+1,"qdagA");
    qdagB=create_2d_double_array(Nr,Ns[1]+1,"qdagB");
    
    // Here is the for loop for setting the propagator initial conditions to 1.0
    for(i=0;i<Nr;i++){
            qA[i][0]=1.0;
            qC[i][0]=1.0;
            //cout<<"A: "<<qintA[i][j]<<endl;
    }

    // Here we solve the diffusion equation for the forwards propagators
    solvediffyQ(qA,w[0],ds,Ns[0],dr);
    
    for(i=0;i<Nr;i++){
        qB[i][0]=qA[i][Ns[0]];
    }
    solvediffyQ(qB,w[1],ds,Ns[1],dr);
    
    //Complementary propagator
    for(i=0;i<Nr;i++){
        qdagB[i][0]=1.0;
    }
    solvediffyQ(qdagB,w[1],ds,Ns[1],dr);
    
    for(i=0;i<Nr;i++){
        qdagA[i][0]=qdagB[i][Ns[1]];
    }
    solvediffyQ(qdagA,w[0],ds,Ns[0],dr);
   
    //Homogeneous propagator
    solvediffyQ(qC,w[2],ds,Ns[2],dr);
    

    
    // Here we get the single chain partition functions Q_AB+Q_C
    Q=q_partition(qdagB,qC,dr,Ns,mu,volume);
        
    cout<<"Q: "<< Q<<endl;
    
    // Here we do the concentration calculation by integration over box and chain
    phi_calc(phi,qA,qdagA,qB,qdagB,qC,Ns,mu,ds);

    //calculation of average concentrations over entire computation box
    phi_total(phi[0],phi[1],phi[2],dr,volume);
    
    
    
    //clearing the memory
    destroy_2d_double_array(qA);
    destroy_2d_double_array(qB);
    destroy_2d_double_array(qC);
    destroy_2d_double_array(qdagA);
    destroy_2d_double_array(qdagB);
    
    return Q;
    
}