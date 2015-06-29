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
    
    
    //Intermediate steps for ADI
    qintA=create_1d_double_array(Nr,"qintA");
    qintB=create_1d_double_array(Nr,"qintB");
    qintC=create_1d_double_array(Nr,"qintC");
    
    
    // Here is the for loop for setting the propagator initial conditions to 1.0
    for(i=0;i<Nr;i++){
            qintA[i]=1.0;
            qA[i][0]=1.0;
           
            qintB[i]=1.0;
            qB[i][0]=1.0;
            
            qintC[i]=1.0;
            qC[i][0]=1.0;
            //cout<<"A: "<<qintA[i][j]<<endl;
    }
    
    
    
    // Here we solve the diffusion equation for the forwards propagators
    solvediffyQ(qA,w[0],qintA,ds,Ns[0],dr);
    solvediffyQ(qB,w[1],qintB,ds,Ns[1],dr);
    solvediffyQ(qC,w[2],qintC,ds,Ns[2],dr);
    
    // The result from the above calculation becomes qdags initial cond
    for(i=0;i<Nr;i++){
                qintA[i]=qB[i][Ns[1]];
                qdagA[i][0]=qB[i][Ns[1]];
                qintB[i]=qA[i][Ns[0]];
                qdagB[i][0]=qA[i][Ns[0]];
                //std::cout<<qintA[i][j][l]<< "----"<<qintB[i][j][l] <<std::endl;
    }
    
    // Here we will solve the diffusion equation for the complementory qs
    solvediffyQ(qdagA,w[0],qintA,ds,Ns[0],dr);
    solvediffyQ(qdagB,w[1],qintB,ds,Ns[1],dr);
    
    // Here we are get the single chain partition functions Q_AB+Q_C
    Q=q_partition(qdagB,qC,dr,Ns,mu);
    

    // Normalizing with respect to box volume
    Q/=volume;
    
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
    destroy_1d_double_array(qintA);
    destroy_1d_double_array(qintB);
    destroy_1d_double_array(qintC);
    
    return Q;
    
};