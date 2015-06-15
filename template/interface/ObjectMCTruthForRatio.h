#ifndef ObjectMCTruthForRatio_H
#define ObjectMCTruthForRatio_H 

#include "math.h"
#include "format_mini.h"
#include "DR.h"
#include <iostream>

/*
   - MC Truth table : 
    -1: unknown 
    1 : W 
    2 : Z 
    3 : from b or B meson
    4 : from c or D meson/bayron
    5 : matched to a parton (light quark)
    15: matched to a gluon 
    nx6: photon
        0(0x6) : prompt photon
        6(1x6) : photon from decay in flight
        12(2x6): photon from ISR 
        18(3x6): photon from FSR 
    7 : tau 
    8 : muon
    9 : electron

 */

int ObjectMCTruthForRatio(int ObjectIdx,GenInfoBranches GenInfo,PhotonInfoBranches ObjectInfo){

        int genidx_=-1;
        float dRmin= 10000;
        for(int g=0;g<GenInfo.Size;g++){
            double dr_ = DR(ObjectInfo.Eta[ObjectIdx],ObjectInfo.Phi[ObjectIdx],GenInfo.Eta[g],GenInfo.Phi[g]);
            double dpT_ = fabs(ObjectInfo.Pt[ObjectIdx]-GenInfo.Pt[g])/GenInfo.Pt[g];
            if(!((GenInfo.Status[g]==1&&GenInfo.PdgID[g]==22)||GenInfo.Status[g]==3)) continue;
            if(fabs(GenInfo.PdgID[g])==23 || fabs(GenInfo.PdgID[g])==24) continue;
            if( dpT_ <3 && dr_ < 0.5 && dr_ < dRmin ){   
                dRmin = dr_;
                genidx_ = g;  
            }    
        }    

        // No matched
        if (genidx_==-1){
            return -1;
            //std::cout<<"ObjectInfo.GenPt : "<<ObjectInfo.GenPt[ObjectIdx]<<" ObjectInfo.Pt : "<<
            //    ObjectInfo.Pt[ObjectIdx]<<" ObjectInfo.Eta : "<<ObjectInfo.Eta[ObjectIdx]<<std::endl;
        // Photon
        }else if(fabs(GenInfo.PdgID[genidx_])==22) {
            if(GenInfo.PhotonFlag[genidx_]==-1){
                return -1;
            }else{
                return GenInfo.PhotonFlag[genidx_]*6.;
            }   
        // b or B meson/baryon
        }else if(fabs(GenInfo.PdgID[genidx_])==5 ||  
                (abs(GenInfo.PdgID[genidx_])%1000)/100 == 5 ||  
                (abs(GenInfo.PdgID[genidx_])%10000)/1000 == 5){ 
            return 3;
        // c or C meson/baryon
        }else if(fabs(GenInfo.PdgID[genidx_])==4 ||  
                (abs(GenInfo.PdgID[genidx_])%1000)/100 == 4 ||  
                (abs(GenInfo.PdgID[genidx_])%10000)/1000 == 4){
            return 4;
        // light quark or gluon
        }else if(fabs(GenInfo.PdgID[genidx_])==21 ||
                (fabs(GenInfo.PdgID[genidx_])>=1&&fabs(GenInfo.PdgID[genidx_])<4) ){
            if(fabs(GenInfo.PdgID[genidx_])==21)
                return 15;
            return 5;
        // Z boson
        }else if(fabs(GenInfo.PdgID[genidx_])==23 ) {
            return 2;
        // W boson
        }else if(fabs(GenInfo.PdgID[genidx_])==24 ) {
            return 1;
        // tau
        }else if(fabs(GenInfo.PdgID[genidx_])==15 ) {
            return 7;
        // muon
        }else if(fabs(GenInfo.PdgID[genidx_])==13 ) {
            return 8;
        // electron
        }else if(fabs(GenInfo.PdgID[genidx_])==11 ) {
            return 9;
        // no considered yet
        }else{
            return -2;
            //std::cout<<"GenInfo.PdgID : "<<fabs(GenInfo.PdgID[genidx_]) <<std::endl;
        }
}

int ObjectMCTruthForRatio(int ObjectIdx,GenInfoBranches GenInfo,LepInfoBranches ObjectInfo){

        int genidx_=-1;
        float dRmin= 10000;
        for(int g=0;g<GenInfo.Size;g++){
            double dr_ = DR(ObjectInfo.Eta[ObjectIdx],ObjectInfo.Phi[ObjectIdx],GenInfo.Eta[g],GenInfo.Phi[g]);
            double dpT_ = fabs(ObjectInfo.Pt[ObjectIdx]-GenInfo.Pt[g])/GenInfo.Pt[g];
            if(!((GenInfo.Status[g]==1&&GenInfo.PdgID[g]==22)||GenInfo.Status[g]==3)) continue;
            if(fabs(GenInfo.PdgID[g])==23 || fabs(GenInfo.PdgID[g])==24) continue;
            if( dpT_ <3 && dr_ < 0.5 && dr_ < dRmin ){   
                dRmin = dr_;
                genidx_ = g;  
            }    
        }    

        // No matched
        if (genidx_==-1){
            return -1;
            //std::cout<<"ObjectInfo.GenPt : "<<ObjectInfo.GenPt[ObjectIdx]<<" ObjectInfo.Pt : "<<
            //    ObjectInfo.Pt[ObjectIdx]<<" ObjectInfo.Eta : "<<ObjectInfo.Eta[ObjectIdx]<<std::endl;
        // Photon
        }else if(fabs(GenInfo.PdgID[genidx_])==22) {
            if(GenInfo.PhotonFlag[genidx_]==-1){
                return -1;
            }else{
                return GenInfo.PhotonFlag[genidx_]*6.;
            }   
        // b or B meson/baryon
        }else if(fabs(GenInfo.PdgID[genidx_])==5 ||  
                (abs(GenInfo.PdgID[genidx_])%1000)/100 == 5 ||  
                (abs(GenInfo.PdgID[genidx_])%10000)/1000 == 5){ 
            return 3;
        // c or C meson/baryon
        }else if(fabs(GenInfo.PdgID[genidx_])==4 ||  
                (abs(GenInfo.PdgID[genidx_])%1000)/100 == 4 ||  
                (abs(GenInfo.PdgID[genidx_])%10000)/1000 == 4){
            return 4;
        // light quark or gluon
        }else if(fabs(GenInfo.PdgID[genidx_])==21 ||
                (fabs(GenInfo.PdgID[genidx_])>=1&&fabs(GenInfo.PdgID[genidx_])<4) ){
            if(fabs(GenInfo.PdgID[genidx_])==21)
                return 15;
            return 5;
        // Z boson
        }else if(fabs(GenInfo.PdgID[genidx_])==23 ) {
            return 2;
        // W boson
        }else if(fabs(GenInfo.PdgID[genidx_])==24 ) {
            return 1;
        // tau
        }else if(fabs(GenInfo.PdgID[genidx_])==15 ) {
            return 7;
        // muon
        }else if(fabs(GenInfo.PdgID[genidx_])==13 ) {
            return 8;
        // electron
        }else if(fabs(GenInfo.PdgID[genidx_])==11 ) {
            return 9;
        // no considered yet
        }else{
            return -2;
            //std::cout<<"GenInfo.PdgID : "<<fabs(GenInfo.PdgID[genidx_]) <<std::endl;
        }
}
int ObjectMCTruthForRatio(int ObjectIdx,GenInfoBranches GenInfo,JetInfoBranches ObjectInfo){

        int genidx_=-1;
        float dRmin= 10000;
        for(int g=0;g<GenInfo.Size;g++){
            double dr_ = DR(ObjectInfo.Eta[ObjectIdx],ObjectInfo.Phi[ObjectIdx],GenInfo.Eta[g],GenInfo.Phi[g]);
            double dpT_ = fabs(ObjectInfo.Pt[ObjectIdx]-GenInfo.Pt[g])/GenInfo.Pt[g];
            if(!((GenInfo.Status[g]==1&&GenInfo.PdgID[g]==22)||GenInfo.Status[g]==3)) continue;
            if(fabs(GenInfo.PdgID[g])==23 || fabs(GenInfo.PdgID[g])==24) continue;
            if( dpT_ <3 && dr_ < 0.5 && dr_ < dRmin ){   
                dRmin = dr_;
                genidx_ = g;  
            }    
        }    

        // No matched
        if (genidx_==-1){
            genidx_=-1;
            dRmin= 10000;
            for(int g=0;g<GenInfo.Size;g++){
                double dr_ = DR(ObjectInfo.Eta[ObjectIdx],ObjectInfo.Phi[ObjectIdx],GenInfo.Eta[g],GenInfo.Phi[g]);
                double dpT_ = fabs(ObjectInfo.Pt[ObjectIdx]-GenInfo.Pt[g])/GenInfo.Pt[g];
                if( dpT_ <3 && dr_ < 0.5 && dr_ < dRmin ){   
                    dRmin = dr_;
                    genidx_ = g;  
                }    
            }    
            //std::cout<<"ObjectInfo.GenPt : "<<GenInfo.Pt[genidx_]<<" ObjectInfo.Pt : "<<
            //    ObjectInfo.Pt[ObjectIdx]<<" ObjectInfo.Eta : "<<ObjectInfo.Eta[ObjectIdx]<<std::endl;
            // ISR/FSR Jet
            if(abs(GenInfo.PdgID[genidx_])==2212||(abs(GenInfo.PdgID[genidx_])<4&&abs(GenInfo.PdgID[genidx_])>0)||
                    GenInfo.PdgID[genidx_]==21){
                if(fabs(GenInfo.PdgID[genidx_])==21)
                    return 15;
                return 5;
            }else{
                return -1;
            }
        // Photon
        }else if(fabs(GenInfo.PdgID[genidx_])==22) {
            if(GenInfo.PhotonFlag[genidx_]==-1){
                return -1;
            }else{
                return GenInfo.PhotonFlag[genidx_]*6.;
            }   
        // b or B meson/baryon
        }else if(fabs(GenInfo.PdgID[genidx_])==5 ||  
                (abs(GenInfo.PdgID[genidx_])%1000)/100 == 5 ||  
                (abs(GenInfo.PdgID[genidx_])%10000)/1000 == 5){ 
            return 3;
        // c or C meson/baryon
        }else if(fabs(GenInfo.PdgID[genidx_])==4 ||  
                (abs(GenInfo.PdgID[genidx_])%1000)/100 == 4 ||  
                (abs(GenInfo.PdgID[genidx_])%10000)/1000 == 4){
            return 4;
        // light quark or gluon
        }else if(fabs(GenInfo.PdgID[genidx_])==21 ||
                (fabs(GenInfo.PdgID[genidx_])>=1&&fabs(GenInfo.PdgID[genidx_])<4) ){
            if(fabs(GenInfo.PdgID[genidx_])==21)
                return 15;
            return 5;
        // Z boson
        }else if(fabs(GenInfo.PdgID[genidx_])==23 ) {
            return 2;
        // W boson
        }else if(fabs(GenInfo.PdgID[genidx_])==24 ) {
            return 1;
        // tau
        }else if(fabs(GenInfo.PdgID[genidx_])==15 ) {
            return 7;
        // muon
        }else if(fabs(GenInfo.PdgID[genidx_])==13 ) {
            return 8;
        // electron
        }else if(fabs(GenInfo.PdgID[genidx_])==11 ) {
            return 9;
        // no considered yet
        }else{
            return -2;
            //std::cout<<"GenInfo.PdgID : "<<fabs(GenInfo.PdgID[genidx_]) <<std::endl;
        }
}

#endif
