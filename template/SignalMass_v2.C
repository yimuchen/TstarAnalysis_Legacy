#include "TStyle.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include <fstream>


void GetMass(TTree *MCTemplatesTree, TH1F *hist);

void SignalMass_v2(){

    gROOT->ProcessLine(".L interface/setTDRStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gStyle->SetErrorX(0.5);


    gStyle->SetOptStat(0);
    const int N = 8;
    int Mass[N] = {650,700,750,800,850,900,950,1000};
    char buffer[128];

    TFile *file = new TFile("MCTemplatesTree.root");
    TH1F *hMass[N];

    /*  
        enum EColor { kWhite =0,   kBlack =1,   kGray=920,
        kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
        kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };
     */
    int hColor[N]={
        800-7,
        600-6,
        432+2,
        616+2,
        900+2,
        kViolet,
        kAzure,
        kTeal
    };

    TCanvas *c1 = new TCanvas("c1","",640,640);

    c1->cd();
    TLegend *legend = new TLegend(0.15,0.75+0.02,0.8+0.02,0.88+0.02+0.02);
    for(int i=0;i<N;i++){
        sprintf(buffer,"tree_Tgamma%iGeV",Mass[i]);
        TTree *ttree_ = (TTree*) file->Get(buffer);
        sprintf(buffer,"hSigMC(%i)",Mass[i]);
        hMass[i] = new TH1F(buffer,"",150,0,1500);
        GetMass(ttree_, hMass[i]);
        hMass[i]->Scale(1./hMass[i]->Integral());
        hMass[i]->SetLineColor(hColor[i]);
        sprintf(buffer,"SigMC(%i)",Mass[i]);
        legend->AddEntry(hMass[i],buffer,"l");
        if(i==0){
            hMass[i]->GetYaxis()->SetRangeUser(0,1.3*hMass[i]->GetMaximum());
            hMass[i]->GetXaxis()->SetTitleOffset(0.95);
            hMass[i]->GetXaxis()->SetTitle("M(top+#gamma)");
            hMass[i]->GetYaxis()->SetTitle("Arbitrary Unit");
            hMass[i]->GetYaxis()->SetTitleOffset(1.4);
            hMass[i]->GetXaxis()->SetLabelSize(0.045);;
            hMass[i]->GetYaxis()->SetLabelSize(0.045);;
            hMass[i]->GetYaxis()->SetTitleSize(0.045);;
            hMass[i]->GetXaxis()->SetTitleSize(0.045);;
            hMass[i]->GetXaxis()->SetNdivisions(505);
            hMass[i]->Draw();
        }else{
            hMass[i]->Draw("same");
            if(i==N-1){
                legend->SetBorderSize(0);
                legend->SetFillColor(0);
                legend->SetFillStyle(0);
                legend->SetNColumns(2);
                legend->SetTextSize(0.04);
                legend->SetTextSizePixels(25);
                legend->Draw();
            }
        }

        delete ttree_;
    }

    c1->SaveAs("signalMassPlot.pdf");


}

void GetMass(TTree *MCTemplatesTree, TH1F *hist){

    float weight_template=-1;
    float Mass_template=-1;
    float category_template=-1;

    MCTemplatesTree->SetBranchAddress("weight",&weight_template);
    MCTemplatesTree->SetBranchAddress("Mass",&Mass_template);
    MCTemplatesTree->SetBranchAddress("category",&category_template);

    for(int entry=0;entry<MCTemplatesTree->GetEntries();entry++){
        MCTemplatesTree->GetEntry(entry);
        if ( weight_template == -1 ) continue;
        if ( Mass_template == -1 ) continue;
        if(category_template==10)
            hist->Fill(Mass_template,weight_template);
    }
}
