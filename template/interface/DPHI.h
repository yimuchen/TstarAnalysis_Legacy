#ifndef DPHI_H
#define DPHI_H


const double PI              = 3.1415926535897932384626433832795;

double DPHI(double phi_1,double phi_2){
    double dphi=-10000;
    if(fabs(phi_1-phi_2)>PI){
        dphi = 2*PI-fabs(phi_1-phi_2);
    }else{
        dphi = fabs(phi_1-phi_2);
    }   
    return dphi;
}

#endif
