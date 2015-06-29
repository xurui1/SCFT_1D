double rterm(int i, double dr, double ds){
    
    if (Coord==1){
        return 0;
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
    }
    
    for (i=0;i<(int)Nr-1;i++){
        rupper[i]=-(ds/(2.0*pow(dr,2.0)))-rterm(i,dr,ds);
        rlower[i]=-(ds/(2.0*pow(dr,2.0)))+rterm(i,dr,ds);
    }

    
    rupper[0]=rupper[0]+rlower[0];
    rlower[Nr-2]=rupper[Nr-2]+rlower[Nr-2];

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
    rupper=create_1d_double_array(Nr-1, "rupper");
    rlower=create_1d_double_array(Nr-1, "rlower");
    
    
    for (s=1;s<(int)Ns+1;s++){
    
        Matrix_r(w,dr,ds,rmid,rupper,rlower);
            
        for (i=0;i<Nr;i++){
            gamma=1.0-(ds/(pow(dr,2.0)))-((ds/2.0)*w[i]);
            betaL=ds/(2.0*pow(dr,2.0))-rterm(i,dr,ds);
            betaU=ds/(2.0*pow(dr,2.0))+rterm(i,dr,ds);
            
            if(i==0){
                bvecr[i]=gamma*qint[i]+(betaL+betaU)*qint[i+1];
            }
            else if(i==(int)Nr-1){
                bvecr[i]=gamma*qint[i]+(betaL+betaU)*qint[i-1];
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
            //cout<<"i: "<<i<<" s: "<<s<<" q: "<<q[i][s]<<endl;
            qint[i]=bvecr[i];
            bvecr[i]=0.0;
        }
        

    }
    
    
    //Deallocate memory
    destroy_1d_double_array(bvecr);
    destroy_1d_double_array(rmid);
    destroy_1d_double_array(rlower);
    destroy_1d_double_array(rupper);


}


