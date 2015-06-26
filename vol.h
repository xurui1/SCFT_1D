double dV(int i,double dr){
    //Definition of volume elements
    double d_V;
    
    if (Coord==1){
        d_V=pow(dr,3.0);
    }
    else if (Coord==2){
        d_V=((double)i*dr)*pow(dr,3.0);
    }
    else if (Coord==3){
        d_V=(pow((double)i*dr,2.0))*pow(dr,3.0);
    }
    
    return d_V;
};