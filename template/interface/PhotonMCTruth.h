#ifndef PhotonMCTRUTH_H
#define PhotonMCTRUTH_H 

#include "math.h"
#include "format.h"
#include "DPHI.h"
#include <iostream>

/*
 Photon from Photon, Jet, Lepton  
              2^2,  2^1, 2^0 
 */

int PhotonMCTruth(int PhotonIdx,GenInfoBranches GenInfo,PhotonInfoBranches PhotonInfo,int &PhotonMatchedGenPdgID){

    int PhotonMCTruth_ = 0;
    int genidx_=-1;
    float dRmin= 10000;
    for(int g=0;g<GenInfo.Size;g++){
        float dphi_ = DPHI(PhotonInfo.GenPhi[PhotonIdx],GenInfo.Phi[g]);
        float deta_ = (PhotonInfo.GenEta[PhotonIdx]-GenInfo.Eta[g]);
        if(  
                fabs(PhotonInfo.GenPt[PhotonIdx]-GenInfo.Pt[g])/GenInfo.Pt[g]<3 &&   
                sqrt(dphi_*dphi_ + deta_*deta_ ) < 0.5 &&
                sqrt(dphi_*dphi_ + deta_*deta_ ) < dRmin
          ){   
            dRmin = sqrt(dphi_*dphi_ + deta_*deta_ );
            genidx_ = g; 
        }    
    }    

    if (genidx_==-1) return 0;
    PhotonMatchedGenPdgID = fabs(GenInfo.PdgID[genidx_]);
    if(fabs(GenInfo.PdgID[genidx_])==22) PhotonMCTruth_ += pow(2.,2); 
    if(fabs(GenInfo.PdgID[genidx_])>=1&&fabs(GenInfo.PdgID[genidx_])<=5) PhotonMCTruth_ += pow(2.,1);
    if(fabs(GenInfo.PdgID[genidx_])>=11&&fabs(GenInfo.PdgID[genidx_])<=16) PhotonMCTruth_ += pow(2.,0);

    if (PhotonMCTruth_==3) std::cout<<"GenInfo.PdgID[genidx_] "<<GenInfo.PdgID[genidx_]<<std::endl;

    return PhotonMCTruth_;
}

#endif
