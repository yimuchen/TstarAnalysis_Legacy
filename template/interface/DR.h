#ifndef DR_H
#define DR_H

const double PI_4dR              = 3.1415926535897932384626433832795;

double DR(double eta_1, double phi_1, double eta_2, double phi_2){
    double dphi=-10000;
    double deta=eta_1-eta_2;
    if(fabs(phi_1-phi_2)>PI_4dR){
        dphi = 2*PI_4dR-fabs(phi_1-phi_2);
    }else{
        dphi = fabs(phi_1-phi_2);
    }   

    return sqrt(dphi*dphi+deta*deta);
}

#endif
