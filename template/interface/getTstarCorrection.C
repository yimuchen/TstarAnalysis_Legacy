#include <iostream>
#include <string>
#include "MCTruthChannel.h"
#include "format.h"
#include "TChain.h"
#include "samples_tr_1lepton4jets.h"
#include <algorithm>
#include <stdio.h>

void getTstarCorrection(){

    int isample_start = Tgamma500GeV;
    int isample_end = TGluonTgamma1500GeV;

    char buffer[256];
    string PATH_="/wk1/cmsdata/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3/2012MC/";
    FILE *pFile;
    pFile = fopen("CorrectionOnTstar.h","w");
    fprintf(pFile,"#ifndef CorrectionOnTstar_H\n");
    fprintf(pFile,"#define CorrectionOnTstar_H\n");
    fprintf(pFile,"float CorrectionOnTstar(int index, int MCTruthChannel_){\n");
    fprintf(pFile,"float MCTruthChannelsDB[%i][4] ={ // [sample][index, hardronic, semilepton, dilepton]\n",
            isample_end-isample_start+1);
    float StandardBR[3] = {4./9, 4./9, 1./9};
    for(int isample_=isample_start;isample_<=isample_end;isample_++){
        string filenameTemp(SAMPLE[isample_].filename);
        filenameTemp.erase(filenameTemp.begin(),filenameTemp.begin()+27);
        filenameTemp.erase(filenameTemp.find('.'),filenameTemp.length());
        sprintf(buffer,"%s%s/*root",PATH_.c_str(),filenameTemp.c_str());
        TChain *root = new TChain("bprimeKit/root");
        root->Add(buffer);

        GenInfoBranches GenInfo;
        GenInfo.Register(root);

        int MCTruthChannels[3] = {0,0,0};
        for(int ievt=0;ievt<root->GetEntries();ievt++){
            root->GetEntry(ievt);
            int MCTruthChannel_ = MCTruthChannel(GenInfo);
            if(MCTruthChannel_!=-1){
                MCTruthChannels[MCTruthChannel_]++;
            }else{
                std::cout<<"[WARNING] NO GEN INFO MATCHED!"<<std::endl;
            }
        }
        std::cout<<filenameTemp <<" ( "<<MCTruthChannels[0]<<" , "<<MCTruthChannels[1]<<" , "<<MCTruthChannels[2]<<" )" <<  std::endl;
        fprintf(pFile,"//%s(%i,%i,%i)\n",filenameTemp.c_str(),MCTruthChannels[0],MCTruthChannels[1],MCTruthChannels[2]);
        if(isample_!=isample_end){
            fprintf(pFile,"{%1.0f, %1.5f, %1.5f, %1.5f},\n",
                    (float) isample_,
                    StandardBR[0]/(MCTruthChannels[0]/(float)root->GetEntries()),
                    StandardBR[1]/(MCTruthChannels[1]/(float)root->GetEntries()),
                    StandardBR[2]/(MCTruthChannels[2]/(float)root->GetEntries()));
        }else{
            fprintf(pFile,"{%1.0f, %1.5f, %1.5f, %1.5f}\n",
                    (float) isample_,
                    StandardBR[0]/(MCTruthChannels[0]/(float)root->GetEntries()),
                    StandardBR[1]/(MCTruthChannels[1]/(float)root->GetEntries()),
                    StandardBR[2]/(MCTruthChannels[2]/(float)root->GetEntries()));
        }
    }
    fprintf(pFile,"};   // MCTruthChannelsDB\n");
    fprintf(pFile,"for(int iscan=0;iscan<%i;iscan++){\n",isample_end-isample_start+1);
    fprintf(pFile,"if(index == (int)MCTruthChannelsDB[iscan][0]){\n");
    fprintf(pFile,"return MCTruthChannelsDB[iscan][MCTruthChannel_+1];\n");
    fprintf(pFile,"}// if\n");
    fprintf(pFile,"}// for\n");
    fprintf(pFile,"return -1;\n");
    fprintf(pFile,"}// CorrectionOnTstar\n");
    fprintf(pFile,"#endif\n");
    fclose(pFile);
}

