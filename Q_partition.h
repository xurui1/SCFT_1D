
double q_partition(double **qdagB,double **qC, double dr, int *Ns, double *mu){
    
    double Q,Q_AB,Q_C;
    int i;
 
    Q=0.0;
    Q_AB=0.0;
    Q_C=0.0;
    
    //ends
    Q_AB+=0.5*qdagB[0][Ns[1]]*dV(0,dr);
    Q_C+=0.5*qC[0][Ns[2]]*dV(0,dr);
    
    Q_AB+=0.5*qdagB[(int)Nr-1][Ns[1]]*dV(Nr-1,dr);
    Q_C+=0.5*qC[(int)Nr-1][Ns[2]]*dV(Nr-1,dr);
    
    //Middle
    for(i=1;i<(int)Nr-1;i++){
            Q_AB+=qdagB[i][Ns[1]]*dV(i,dr);
            Q_C+=qC[i][Ns[2]]*dV(i,dr);
    }
    
    Q_AB=exp(mu[0])*Q_AB;
    Q_C=(exp(mu[1]*kappa)*Q_C)/kappa;
    cout<<"Q_AB: "<<Q_AB<<" Q_C: "<<Q_C<<endl;
    //I'm adding the two single chain partition functions together for the return function
    Q=Q_AB+Q_C;
    
    
    return Q;
};