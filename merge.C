#include "template/interface/samples_tr_1lepton4jets.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

// Normal/MCTemplatesTree.root
//   KEY: TTree tree_Tgamma1500GeV;1    tree_Tgamma1500GeV
// tree_[SAMPLE[samples_order_size].tag][RunStatusNames[RunStatus]]
// KEY: TTree    ResultTree;1    ResultTree


char buffer[256];
void merge(){

    TFile *files[(const int) (RunStatusSize)];
    TTree *trees[(const int)(samples_order_size*RunStatusSize)];
    TTree *ResultTree;
    TChain *ResultTree_tmp = new TChain("ResultTree");


    TFile *output = new TFile("SumNtuples.root","recreate");
    for(int iunc=0;iunc<RunStatusSize;iunc++){
        sprintf(buffer,"%s/MCTemplatesTree.root",RunStatusNames[iunc].c_str());
        if(iunc==Normal)
            sprintf(buffer,"Normal/MCTemplatesTree.root");
        printf("[Loading] %s\n",buffer);
        files[iunc] = new TFile(buffer);
        ResultTree_tmp->Add(buffer);

        for(int isample = 0; isample < samples_order_size; isample++){
            sprintf(buffer,"tree_%s%s",SAMPLE[isample].tag, RunStatusNames[iunc].c_str());
            TTree *tree_tmp = (TTree*) files[iunc]->Get(buffer);
            if(tree_tmp){
                output->cd();

                trees[iunc*samples_order_size + isample] = tree_tmp->CloneTree(0);
                for(int ientry = 0; ientry < tree_tmp->GetEntries(); ientry++){
                    tree_tmp->GetEntry(ientry);
                    trees[iunc*samples_order_size + isample]->Fill();
                }
                trees[iunc*samples_order_size + isample]->Write();
            }
            delete tree_tmp;

        }

    }
    if(ResultTree_tmp){
        ResultTree = ResultTree_tmp->CloneTree(0);
        output->cd();
        for(int ientry = 0; ientry < ResultTree_tmp->GetEntries(); ientry++){
            ResultTree_tmp->GetEntry(ientry);
            ResultTree->Fill();
        }
    }
    delete ResultTree_tmp;
    ResultTree->Write();
    output->Close();
}
