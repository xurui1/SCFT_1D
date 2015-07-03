void omega(double **w){
    
    int i,x;
    
    //This is for a sinusoidal initial omega field. Adjust 'initial' in the parameters header.
    //Currently will give a lamellar structure in z-direction
    if (initial==0){
        for(i=0;i<Nr;i++){
            w[0][i]=-5.0*cos(2.0*Pi*(double)i/(double)Nr)+(double)((rand() %100)/50.0)-1.0;
            w[1][i]=5.0*cos(2.0*Pi*(double)i/(double)Nr)+(double)((rand() %100)/50.0)-1.0;
            w[2][i]=10.0+(double)((rand() %100)/50.0-1.0);
        }
    }
    
    
    //This is for a bilayer conformation
    else if (initial==1){
        ifstream Init1;
        Init1.open("bilayer_M50_N50.dat");
        for (i=0;i<Nr;i++){
                Init1 >> x >> w[0][i]>>w[1][i]>>w[2][i];
        }
        Init1.close();
    }
    
    for (i=0;i<Nr;i++){
                std::cout << w[0][i]<<" "<<w[1][i]<<" "<<w[2][i]<<endl;
    }
    
}