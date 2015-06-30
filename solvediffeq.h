double rterm(int i, double dr, double ds){
    
    if (Coord==1){
        return 0.0;
    }
    else{
        return ((double)Coord-1.0)*(ds/((dr*((double)i)+r_0)*4.0*dr));
    }
}




//Build LHS matrix
void Matrix_r(double *w, double dr, double ds, double *rmid, double *rupper, double *rlower){
    int i;
    
    for (i=0;i<Nr;i++){
        rmid[i]=1.0+(ds/(pow(dr,2.0)))+(ds/2.0)*w[i];
        rupper[i]=-(ds/(2.0*pow(dr,2.0)))-rterm(i-1,dr,ds);
        rlower[i]=-(ds/(2.0*pow(dr,2.0)))+rterm(i+1,dr,ds);
    }
    
    rupper[0]=rupper[0]+rlower[0];
    rlower[Nr-1]=rupper[Nr-1]+rlower[Nr-1];
    rupper[Nr-1]=0.0;
    rlower[0]=0.0;

}



//Apply Crank-Nicolson method to solve the modified diffusion equation
void solvediffyQ(double **q, double *w, double *qint, double ds, int Ns, double dr){

    int            i,s;  // some counters
    double gamma, betaL, betaU;
    double *bvecr;
    double *rmid;
    double *rupper;
    double *rlower;
    
    //Allocate memory for `matrix' operations
    bvecr=create_1d_double_array(Nr, "bvecr");
    rmid=create_1d_double_array(Nr, "rmid");
    rupper=create_1d_double_array(Nr, "rupper");
    rlower=create_1d_double_array(Nr, "rlower");
    
    
    for (s=1;s<(int)Ns+1;s++){
    
        Matrix_r(w,dr,ds,rmid,rupper,rlower);
            
        for (i=0;i<Nr;i++){
            gamma=1.0-(ds/(pow(dr,2.0)))-((ds/2.0)*w[i]);
            betaL=ds/(2.0*pow(dr,2.0))-rterm(i-1,dr,ds);
            betaU=ds/(2.0*pow(dr,2.0))+rterm(i+1,dr,ds);
            //cout<<"i: "<<i<<" gamma: "<<gamma<<" betaL: "<<betaL<<" betaU: "<<betaU<<endl;
            
            if(i==0){
                bvecr[i]=gamma*qint[0]+(betaL+betaU)*qint[1];
            }
            else if(i==(int)Nr-1){
                bvecr[i]=gamma*qint[Nr-1]+(betaL+betaU)*qint[Nr-2];
            }
            else{
                bvecr[i]=gamma*qint[i]+betaL*qint[i-1]+betaU*qint[i+1];
            }
        }
        
        //Apply TDMA
        TDMA(bvecr,Nr,rlower,rmid,rupper);
            
        //Now we have our solution for all i for s, from s-1. Full step completed.
        for (i=0;i<Nr;i++){
            q[i][s]=bvecr[i];
            qint[i]=bvecr[i];
            bvecr[i]=0.0;
            if (fabs(q[i][s])>1e4){
                cout<<i<<" "<<s<<" propagator problem: "<<q[i][s]<<endl;
                exit(EXIT_FAILURE);
            }
            
        }
        

    }
    
    
    //Deallocate memory
    destroy_1d_double_array(bvecr);
    destroy_1d_double_array(rmid);
    destroy_1d_double_array(rlower);
    destroy_1d_double_array(rupper);


}


