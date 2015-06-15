#ifndef RecoNoCutPhotonSelection_H
#define RecoNoCutPhotonSelection_H

#include "format.h"
#include "DPHI.h"
#include "TLorentzVector.h"

void RecoNoCutPhotonSelection(PhotonInfoBranches PhotonInfo, 
        int &NPhotons, int *P_Index){

    for(int nl=0;nl<PhotonInfo.Size;nl++){
        // Loose Photon
        if(     
                (fabs(PhotonInfo.Pt[nl])>30) &&
                (// Barrel
                  fabs(PhotonInfo.Eta[nl])<=1.4442 
                 ||
                 // Endcap
                  (fabs(PhotonInfo.Eta[nl])<2.1)&&(fabs(PhotonInfo.Eta[nl])>=1.566) 
                 )
          ){
            P_Index[NPhotons] = nl;
            NPhotons++;
        }
    }   

}

#endif
