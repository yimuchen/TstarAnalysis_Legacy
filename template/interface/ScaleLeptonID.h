#ifndef ScaleLeptonID_H
#define ScaleLeptonID_H

#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include <iostream>

bool debug = false;
TFile *file = new TFile("interface/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root");
TGraphAsymmErrors *muonid = (TGraphAsymmErrors*)file->Get("DATA_over_MC_combRelIsoPF04dBeta<012_Tight_eta_pt20-500_2012ABCD");

float ScaleLeptonID(int leptonType,float pT, float eta, int MODE){
    float scale_ = -1;
    /*
       electron ID :
            CutBased
            https://twiki.cern.ch/twiki/bin/view/Main/EGammaScaleFactors2012#2012_8_TeV_data_53X

            MVA
            https://twiki.cern.ch/twiki/bin/view/CMS/KoPFAElectronTagAndProbe
       muon ID : 
https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs from Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root
     */
    if(leptonType == 11){
        // [3] [3] = [eta] [pT]
        // MVA > 0.9
        float scale[3][3] = { 
            {0.939, 0.950, 0.957},
            {0.920, 0.949, 0.959},
            {0.907, 0.937, 0.954}
        };  
        float scaleErr[3][3] = { 
            {0.003, 0.001, 0.001},
            {0.002, 0.002, 0.003},
            {0.005, 0.008, 0.011}
        };  

        float etaRegion[3][2] = { 
            {0.0,0.8},
            {0.8,1.478},
            {1.478,2.5}
        };  
        float pTRegion[3][2] = { 
            {30,40},
            {40,50},
            {50,100000000}
        };  

        for(int i_eta = 0;i_eta<3;i_eta++)
        for(int i_pT = 0;i_pT<3;i_pT++){
            if(pT>=pTRegion[i_pT][0] && pT<pTRegion[i_pT][1])
                if(fabs(eta)>=etaRegion[i_eta][0] && fabs(eta)<etaRegion[i_eta][1]){
                    scale_ = scale[i_eta][i_pT] + (MODE - 1.)*scaleErr[i_eta][i_pT];
                    break;
                }   
        }   

    }else if(leptonType == 13){

        const int N = muonid->GetN();
        float etaRegion[N][2]; 
        float scale[N];
        float scaleErr[N];
        for(int i=0;i<N;i++){
            double Yscale = 0;
            double Xeta = 0;
            muonid->GetPoint(i,Xeta,Yscale);
            scale[i] = Yscale;
            etaRegion[i][0] = Xeta-muonid->GetErrorXlow(i);
            etaRegion[i][1] = Xeta+muonid->GetErrorXhigh(i);
            scaleErr[i] = muonid->GetErrorYlow(i);
            if(muonid->GetErrorYhigh(i)>scaleErr[i])
                scaleErr[i] = muonid->GetErrorYhigh(i);
            if(debug)
                std::cout<<i<<"th ( Xeta, Yscale) = ( "<<etaRegion[i][0] <<" ~ "<< 
                    etaRegion[i][1] <<" , "<<scale[i]<<" ) "<<std::endl;
        }   
        for(int i_eta = 0;i_eta<N;i_eta++)
            if(eta>=etaRegion[i_eta][0] && eta<etaRegion[i_eta][1]){
                scale_ = scale[i_eta] + (MODE - 1.)*scaleErr[i_eta];
                break;
            }
    }

    return scale_;
}

#endif
