#include "proj_anyregions_shape_loose_norunCR.cc"
#include "interface/ConstantNumbers.h"
#include "interface/samples_tr_1lepton4jets.h"
#include <iostream>
#include "TChain.h"
#include "pdf_shifttable_log__STInSR.h"
#include "pdf_shifttable_log__ST.h"

float UncOnLumi = 2.6; //https://hypernews.cern.ch/HyperNews/CMS/get/B2G-12-014/70.html
float UncOnTrig = 1.0; //https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Trigger_Efficiencies
float UncOnISR[2] = {2.72, 7.27};   // [Control, Signal] from T2TopGamma/Analysis/T2TopPlusPhoton/ISR_FSR
float UncOnFSR[2] = {1.53, 5.65};   // [Control, Signal] from T2TopGamma/Analysis/T2TopPlusPhoton/ISR_FSR
const int NSigAndBkg_ = samples_order_size;
float TotalUncs_[NSigAndBkg_];
int SignalSamples[2] = { Tgamma700GeV, Tgamma1200GeV};
//int SignalSamples[2] = { Tgamma500GeV, Tgamma1500GeV};
int SignalRegionIndex = _MtstarChi2_2r1l4j;
const int Nuncs = 11+2+1+1;
string Uncs_names[Nuncs] = {
    "Lumi",
    "Trigger",
    "JES",
    "JER",
    "LepID",
    "PhoID",
    "PU",
    "Xsec",
    "Scale",
    "Matching",
    "ISR",
    "FSR",
    "TopPt",
    "PDF",
    "Statistics"
};
string Uncs_Input_names[Nuncs] = {
    "2.6\\\%",
    "Top Group",
    "JetMET Group",
    "JetMET Group",
    "POG ID",
    "POG ID",
    "$\\sigma_{\\rm MiniBias}$",
    "CMS or theory",
    "$\\ttbar/DY+$jets sample ",
    "$\\ttbar/DY+$jets sample ",
    "ISR",
    "FSR",
    "Top $\\pt$ reweighting",
    "CTEQ61",
    "MC"
};
enum Uncs_numb{
    _Lumi,
    _Trigger,
    _JES,
    _JER,
    _LepID,
    _PhoID,
    _PU,
    _Xsec,
    _Scale,
    _Matching,
    _ISR,
    _FSR,
    _TopPt,
    _PDF,
    _statistics,
    _Uncs_numb_size
};

void JES_JER_PU(TChain *ResultTree, int RegionIndex, float unc[], int UNC_MODE);
void Stat(TChain *ResultTree, int RegionIndex, float unc[]);
void Xsec(TChain *ResultTree, int RegionIndex, float unc[]);
void Scale_Match(TChain *ResultTree, int RegionIndex, float unc[], int UNC_MODE);
void Lep_Pho_ID(TChain *ResultTree, int RegionIndexB4, int RegionIndexaf, float unc[]);
void Display(TChain *ResultTree, int _RegionNoLSPS, int _RegionNoPS, int _Region,float &syst_bkg);
void TableForSFbySF(TChain *ResultTree);
void TableForYield(TChain *ResultTree, int Region_, float syst_bkg);
void TrimChar(char *OrignalChar_, char *TrimChar_);

void MakeSystematics(){

    TChain *ResultTree = new TChain("ResultTree");
    ResultTree->Add("SumNtuples.root");

    float syst_bkgInCR = 0.;
    float syst_bkgInSR = 0.;

    printf("\n");
    printf("Control Region\n");
    Display(ResultTree,_NOfVertexafPUTrgeff1p1l4j,_NOfVertexafPUTrgeff1p1l4jafLS,_NOfVertexafPUTrgeff1p1l4jafLSPSTop,syst_bkgInCR);

    printf("\n");
    printf("Signal Region\n");
    Display(ResultTree, _MtstarChi2_2r1l4jnoLSPS, _MtstarChi2_2r1l4jnoPS, _MtstarChi2_2r1l4j,syst_bkgInSR);

    /*
       For Latex
     */
    TableForSFbySF(ResultTree);
    TableForYield(ResultTree, _NOfVertexafPUTrgeff1p1l4jafLSPSTop, syst_bkgInCR);
    TableForYield(ResultTree, _MtstarChi2_2r1l4j, syst_bkgInSR);

    //std::cout<<"(syst_bkgInCR, syst_bkgInSR) = ("<<syst_bkgInCR<<" , "<<syst_bkgInSR<<")"<<std::endl;

    printf("\n\n// [Signal Unc for EXO limit]\n double SignalTotalUnc[%i][2] = {\n", SignalSamples[1]-SignalSamples[0]+1);
    for(int isample=1;isample<SignalSamples[1]-SignalSamples[0]+2;isample++){
        char TrimChar_[128];
        TrimChar(SAMPLE[SignalSamples[0]+isample-1].tag, TrimChar_);
        if(isample!= SignalSamples[1]-SignalSamples[0]+1){
            printf("{%i, %13.1f},\n ",atoi(TrimChar_),sqrt(TotalUncs_[isample])); 
        }else{
            printf("{%i, %13.1f}};\n ",atoi(TrimChar_),sqrt(TotalUncs_[isample])); 
        }
    }
}

void TableForSFbySF(TChain *ResultTree){
    int ResultTree_plot = -1;
    int ResultTree_sampleSize = samples_order_size + 1; // for total MC
    int ResultTree_sampleID[ResultTree_sampleSize];
    double ResultTree_yield[ResultTree_sampleSize];
    double ResultTree_stat_error[ResultTree_sampleSize];
    for(int isample=0;isample<ResultTree_sampleSize;isample++){
        ResultTree_sampleID[isample] = -2;
        ResultTree_yield[isample] = -2;
        ResultTree_stat_error[isample] = -2;
    }    
    int ResultTree_RunMode = -1;
    ResultTree->SetBranchAddress("ResultTree_sampleSize",&ResultTree_sampleSize);
    ResultTree->SetBranchAddress("ResultTree_plot",&ResultTree_plot);
    ResultTree->SetBranchAddress("ResultTree_sampleID",&ResultTree_sampleID[0]);
    ResultTree->SetBranchAddress("ResultTree_yield",&ResultTree_yield[0]);
    ResultTree->SetBranchAddress("ResultTree_stat_error",&ResultTree_stat_error[0]);
    ResultTree->SetBranchAddress("ResultTree_RunMode",&ResultTree_RunMode);
    float Yield_SFbySF[7][2];
    int index_[6] = {
        _NOfVertexb4PU1p1l4j,
        _NOfVertexafPU1p1l4j,
        _NOfVertexafPUTrgeff1p1l4j,
        _NOfVertexafPUTrgeff1p1l4jafLS,
        _NOfVertexafPUTrgeff1p1l4jafLSPS,
        _NOfVertexafPUTrgeff1p1l4jafLSPSTop
    };
    for(int itree=0;itree<ResultTree->GetEntries();itree++){
        ResultTree->GetEntry(itree);
        if(ResultTree_RunMode != Normal) continue;

        for(int i=0;i<6;i++){
            if(ResultTree_plot == index_[i]){
                Yield_SFbySF[i][0] = ResultTree_yield[0];
                Yield_SFbySF[i][1] = ResultTree_stat_error[0];

                // data
                Yield_SFbySF[6][0] = ResultTree_yield[1];
                Yield_SFbySF[6][1] = ResultTree_stat_error[1];
            }
        }
    }
    printf("\n------- Table for SF by SF ----------\n");
    printf("\\begin{tabular}{|c|c|c|c|c|c|c|}\n");
    printf("\\hline\n");
    printf("& \\multicolumn{6}{c|}{Yields} \\\\ \n");
    printf("& no scale factors & + pile-up & + trigger & + lepton Id & + photon Id & + Top $\\pt$  \\\\ \n");
    printf("\\hline\n");
    printf("Total bkg. (MC) ");
    for(int i=0;i<6;i++)
        printf("& %1.0f $\\pm$ %1.0f ",Yield_SFbySF[i][0], Yield_SFbySF[i][1]);
    printf("\\\\ \n");
    printf("\\hline\n");
    printf("data & \\multicolumn{6}{c|}{%1.0f $\\pm$ %1.0f} \\\\ \n",Yield_SFbySF[6][0], Yield_SFbySF[6][1]);
    printf("\\hline\n");
    printf("\\end{tabular} \n");

}
void TableForYield(TChain *ResultTree, int Region_, float syst_bkg){
    int ResultTree_plot = -1;
    int ResultTree_sampleSize = samples_order_size + 1; // for total MC
    int ResultTree_sampleID[ResultTree_sampleSize];
    double ResultTree_yield[ResultTree_sampleSize];
    double ResultTree_stat_error[ResultTree_sampleSize];
    for(int isample=0;isample<ResultTree_sampleSize;isample++){
        ResultTree_sampleID[isample] = -2;
        ResultTree_yield[isample] = -2;
        ResultTree_stat_error[isample] = -2;
    }    
    int ResultTree_RunMode = -1;
    ResultTree->SetBranchAddress("ResultTree_sampleSize",&ResultTree_sampleSize);
    ResultTree->SetBranchAddress("ResultTree_plot",&ResultTree_plot);
    ResultTree->SetBranchAddress("ResultTree_sampleID",&ResultTree_sampleID[0]);
    ResultTree->SetBranchAddress("ResultTree_yield",&ResultTree_yield[0]);
    ResultTree->SetBranchAddress("ResultTree_stat_error",&ResultTree_stat_error[0]);
    ResultTree->SetBranchAddress("ResultTree_RunMode",&ResultTree_RunMode);
    const int Nsamples_ = 41;
    float YieldForTable[Nsamples_][2];
    int index_[Nsamples_] = {
        // Signal
        Tgamma700GeV,
        Tgamma750GeV,
        Tgamma800GeV,
        Tgamma850GeV,
        Tgamma900GeV,
        Tgamma950GeV,
        Tgamma1000GeV,
        Tgamma1100GeV,
        Tgamma1200GeV,
        Tgamma1300GeV,
        Tgamma1400GeV,
        Tgamma1500GeV,
        // Bosons
        DYJetsToLL_M_10To50filter_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
        DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1,
        WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2,
        WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        WWWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
        WGToLNuG_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        WWGJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
        ZGToLLG_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
        GJets_HT_40To100_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V19_v1,
        GJets_HT_100To200_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V19_v1,
        GJets_HT_200To400_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
        GJets_HT_400ToInf_8TeV_madgraph_v3_Summer12_DR53X_PU_S10_START53_V7C_v1,
        DiPhotonJets_8TeV_madgraph_tarball_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
        // ttbar
        TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        // Single Top
        T_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        T_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v2,
        T_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        Tbar_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        Tbar_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        Tbar_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
        // tt+X
        TTGJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
        TTWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
        TTWWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
        TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
        TTH_Inclusive_M_125_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
        // Total bkg.
        -1,
        // data
        Data
    };
    for(int itree=0;itree<ResultTree->GetEntries();itree++){
        ResultTree->GetEntry(itree);
        if(ResultTree_RunMode != Normal) continue;
        if(ResultTree_plot != Region_) continue;

        for(int isample_=0;isample_<Nsamples_;isample_++){
            YieldForTable[isample_][0] = ResultTree_yield[1+index_[isample_]];
            YieldForTable[isample_][1] = ResultTree_stat_error[1+index_[isample_]];
        }

    }
    printf("\n------- Table for Yields in %s ----------\n", PLOTS[Region_].PlotName);
    printf("\\begin{tabular}{|l|l|c|}\n");
    printf("\\hline\n");
    printf("\\multicolumn{2}{|c|}{Process} & Yields \\\\ \n");
    printf("\\hline\n");
    printf("\\multirow{12}{*}{Signal}\n");

    for(int isample_=0;isample_<Nsamples_;isample_++){
        if(index_[isample_]==DYJetsToLL_M_10To50filter_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1){
            printf("\\hline\n");
            printf("\\hline\n");
            printf("\\multirow{15}{*}{Boson(s)}\n");
        }else if (index_[isample_]==TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1){
            printf("\\hline\n");
        }else if (index_[isample_]==T_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1){
            printf("\\hline\n");
            printf("\\multirow{6}{*}{Single top}\n");
        }else if (index_[isample_]==TTGJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1){
            printf("\\hline\n");
            printf("\\multirow{5}{*}{\\ttbar+X}\n");
        }else if (index_[isample_]==-1){
            printf("\\hline\n");
        }else if (index_[isample_]==Data){
            printf("\\hline\n");
            printf("\\hline\n");
        }


        if(index_[isample_]!=-1&&
                index_[isample_]!=TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1&&
                index_[isample_]!=Data){
            if(YieldForTable[isample_][0]!=0){
                printf("& %s & %1.3f $\\pm$ %1.3f \\\\ \n",
                        SAMPLE[index_[isample_]].latex,
                        YieldForTable[isample_][0],
                        YieldForTable[isample_][1]);
            }else{
                printf("& %s & $<$ %1.3f \\\\ \n",
                        SAMPLE[index_[isample_]].latex,
                        YieldForTable[isample_][1]);
            }
        }else{
            if(index_[isample_]==-1){
                printf("\\multicolumn{2}{|c|}{Total bkg. (MC)} & %1.3f $\\pm$ %1.3f (syst.) \\\\ \n",
                        YieldForTable[isample_][0],
                        YieldForTable[isample_][0]*syst_bkg/100.);
                        //YieldForTable[isample_][1]);
            }else{
                printf("\\multicolumn{2}{|c|}{%s} & %1.3f $\\pm$ %1.3f \\\\ \n",
                        SAMPLE[index_[isample_]].latex,
                        YieldForTable[isample_][0],
                        YieldForTable[isample_][1]);
            }
        }
    }
    printf("\\hline\n");
    printf("\\end{tabular}\n");


}

void Display(TChain *ResultTree, int _RegionNoLSPS, int _RegionNoPS, int _Region, float &syst_bkg){
    int IsSR = 1;
    if(_Region!=SignalRegionIndex) IsSR = 0;

    float Uncs[Nuncs][SignalSamples[1]-SignalSamples[0]+2];
    float TotalUncs[SignalSamples[1]-SignalSamples[0]+2];
    for(int idisp = 0;idisp<Nuncs;idisp++)
        for(int isample=0;isample<SignalSamples[1]-SignalSamples[0]+2;isample++){
            Uncs[idisp][isample] = 0.0; 
            TotalUncs[isample] = 0.0; 
        }

    for(int isample = 0;isample<SignalSamples[1]-SignalSamples[0]+2;isample++) {
        Uncs[_Lumi][isample]    = UncOnLumi;
        //Uncs[_Trigger][isample] = UncOnTrig;
    }

    JES_JER_PU(ResultTree, _Region,Uncs[_JES],_JES);
    JES_JER_PU(ResultTree, _Region,Uncs[_JER],_JER);
    JES_JER_PU(ResultTree, _Region,Uncs[_PU],_PU);
    Xsec(ResultTree, _Region, Uncs[_Xsec]);
    Stat(ResultTree, _Region, Uncs[_statistics]);
    Scale_Match(ResultTree, _Region, Uncs[_Scale],_Scale);
    Scale_Match(ResultTree, _Region, Uncs[_Matching],_Matching);
    //Lep_Pho_ID(ResultTree, _RegionNoLSPS, _RegionNoPS, Uncs[_LepID]);
    JES_JER_PU(ResultTree, _Region,Uncs[_LepID],_LepID);
    //Lep_Pho_ID(ResultTree, _RegionNoPS, _Region, Uncs[_PhoID]);
    JES_JER_PU(ResultTree, _Region,Uncs[_PhoID],_PhoID);
    JES_JER_PU(ResultTree, _Region,Uncs[_TopPt],_TopPt);
    JES_JER_PU(ResultTree, _Region,Uncs[_Trigger],_Trigger);

    for(int isample = 1;isample<SignalSamples[1]-SignalSamples[0]+2;isample++) {
        Uncs[_ISR][isample] = UncOnISR[IsSR];
        Uncs[_FSR][isample] = UncOnFSR[IsSR];
    }

    for(int isample = 0;isample<SignalSamples[1]-SignalSamples[0]+2;isample++) {
        if(IsSR){
            if(isample == 0){
                Uncs[_PDF][isample] = fabs(overall_background__STInSR_plus);
                if(Uncs[_PDF][isample]<fabs(overall_background__STInSR_minus)) 
                    Uncs[_PDF][isample] = fabs(overall_background__STInSR_minus);
            }else{
                Uncs[_PDF][isample] = fabs(PDF_unc__STInSR_plus[isample+SignalSamples[0]-1]);
                if(Uncs[_PDF][isample]<fabs(PDF_unc__STInSR_minus[isample+SignalSamples[0]-1]))
                    Uncs[_PDF][isample] = fabs(PDF_unc__STInSR_minus[isample+SignalSamples[0]-1]);
            }
            /*
            printf("[In SR] %ith sample bkg(+%f, -%f), sig(+%f, -%f), Unc(%f)\n",
                    isample,
                    fabs(overall_background__STInSR_plus), fabs(overall_background__STInSR_minus),
                    fabs(PDF_unc__STInSR_plus[isample+SignalSamples[0]-1]),
                    fabs(PDF_unc__STInSR_minus[isample+SignalSamples[0]-1]),
                    Uncs[_PDF][isample]
                    );
            */
        }else{
            if(isample == 0){
                Uncs[_PDF][isample] = fabs(overall_background__ST_plus);
                if(Uncs[_PDF][isample]<fabs(overall_background__ST_minus)) 
                    Uncs[_PDF][isample] = fabs(overall_background__ST_minus);
            }else{
                Uncs[_PDF][isample] = fabs(PDF_unc__ST_plus[isample+SignalSamples[0]-1]);
                if(Uncs[_PDF][isample]<fabs(PDF_unc__ST_minus[isample+SignalSamples[0]-1]))
                    Uncs[_PDF][isample] = fabs(PDF_unc__ST_minus[isample+SignalSamples[0]-1]);
            }
        }
        Uncs[_PDF][isample] *= 100.; // convert to unit (100%)
    }

    for(int idisp = 0;idisp<_Uncs_numb_size;idisp++){
        if(idisp==0) {
            printf("          |%13s|","Total MC"); 
            for(int isample=1;isample<SignalSamples[1]-SignalSamples[0]+2;isample++)
                printf("%13s|",SAMPLE[SignalSamples[0]+isample-1].tag); 
            printf("\n"); 

            for(int iloop = 0;iloop<10;iloop++)
                printf("-");
            printf("|");
            for(int isample=0;isample<SignalSamples[1]-SignalSamples[0]+2;isample++){
                for(int iloop = 0;iloop<13;iloop++)
                    printf("-");
                printf("|");
            }
            printf("\n"); 
        }
        printf("%10s|",Uncs_names[idisp].c_str());
        for(int isample=0;isample<SignalSamples[1]-SignalSamples[0]+2;isample++){
            TotalUncs[isample] += Uncs[idisp][isample]*Uncs[idisp][isample];
            if(IsSR){
                printf("%13.1f|",Uncs[idisp][isample]); 
            }else{
                printf("%13.2f|",Uncs[idisp][isample]); 
            }
        }
        printf("\n"); 

    }
    for(int iloop = 0;iloop<10;iloop++)
        printf("-");
    printf("|");
    for(int isample=0;isample<SignalSamples[1]-SignalSamples[0]+2;isample++){
        for(int iloop = 0;iloop<13;iloop++)
            printf("-");
        printf("|");
    }
    printf("\n"); 
    printf("%10s|","Total");
    for(int isample=0;isample<SignalSamples[1]-SignalSamples[0]+2;isample++){
        if(IsSR){
            printf("%13.1f|",sqrt(TotalUncs[isample])); 
        }else{
            printf("%13.2f|",sqrt(TotalUncs[isample])); 
        }
    }
    printf("\n"); 

    // For Latex
    char TrimChar_[128];
    printf("\n\n for Latex \n\n"); 

    printf("\\begin{tabular}{|l||l|c|");
    for(int isample=0;isample<SignalSamples[1]-SignalSamples[0]+1;isample++)
        printf("c");
    printf("|}\n");
    printf("\\hline\n");
    for(int idisp = 0;idisp<_Uncs_numb_size;idisp++){
        if(idisp==0) {
            printf("\\multirow{2}{*}{Source} & \\multirow{2}{*}{Input} & \\multirow{2}{*}{Total Bkg.} & \\multicolumn{%i}{c|}{ $\\cPqtstar\\overline{\\cPqtstar}$, $M(\\cPqtstar) [\\GeVcc] =$} \\\\ \n", SignalSamples[1]-SignalSamples[0]+1); 
            printf(" &  &  "); 
            for(int isample=1;isample<SignalSamples[1]-SignalSamples[0]+2;isample++){
                TrimChar(SAMPLE[SignalSamples[0]+isample-1].tag, TrimChar_);
                printf("& %13s ",TrimChar_); 
            }
            printf("\\\\ \n"); 
            printf("\\hline\n");

        }
        printf("%10s ",Uncs_names[idisp].c_str());
        printf("& %10s ",Uncs_Input_names[idisp].c_str());

        for(int isample=0;isample<SignalSamples[1]-SignalSamples[0]+2;isample++){
            if(Uncs[idisp][isample]==0.){
                printf("& - "); 
            }else{
                if(IsSR){
                    printf("& %13.1f ",Uncs[idisp][isample]); 
                }else{
                    printf("& %13.2f ",Uncs[idisp][isample]); 
                }
            }
        }
        printf("\\\\ \n"); 

    }
    printf("\\hline\n");
    printf("\\hline\n");
    printf("%10s & ","Sum [\\%]");
    for(int isample=0;isample<SignalSamples[1]-SignalSamples[0]+2;isample++){
        if(IsSR){
            printf("& %13.1f ",sqrt(TotalUncs[isample])); 
        }else{
            printf("& %13.2f ",sqrt(TotalUncs[isample])); 
        }
        TotalUncs_[isample] = TotalUncs[isample];
    }
    printf("\\\\ \n"); 
    printf("\\hline\n");
    printf("\\end{tabular}\n");
    printf("\n"); 
    printf("\n"); 

    syst_bkg = sqrt(TotalUncs[0]);
}

void JES_JER_PU(TChain *ResultTree, int RegionIndex, float unc[], int UNC_MODE){
    int RunModes_[2];
    if(UNC_MODE == _JES){
        RunModes_[0] = UncJESMinus;
        RunModes_[1] = UncJESPlus;
    }else if(UNC_MODE == _JER){
        RunModes_[0] = UncJERMinus;
        RunModes_[1] = UncJERPlus;
    }else if(UNC_MODE == _PU){
        RunModes_[0] = UncPUMinus;
        RunModes_[1] = UncPUPlus;
    }else if(UNC_MODE == _PhoID){
        RunModes_[0] = UncPhoIDMinus;
        RunModes_[1] = UncPhoIDPlus;
    }else if(UNC_MODE == _LepID){
        RunModes_[0] = UncLepIDMinus;
        RunModes_[1] = UncLepIDPlus;
    }else if(UNC_MODE == _Trigger){
        RunModes_[0] = UncTrigMinus;
        RunModes_[1] = UncTrigPlus;
    }else if(UNC_MODE == _TopPt){
        RunModes_[0] = UncTopPtMinus;
        RunModes_[1] = UncTopPtPlus;
    }else{
        std::cout<<"[ERROR] This function is not for UNC_MODE of "<<UNC_MODE<<std::endl;
        exit(0);
    }
    int ResultTree_plot = -1;
    int ResultTree_sampleSize = samples_order_size + 1; // for total MC
    int ResultTree_sampleID[ResultTree_sampleSize];
    double ResultTree_yield[ResultTree_sampleSize];
    double ResultTree_stat_error[ResultTree_sampleSize];
    for(int isample=0;isample<ResultTree_sampleSize;isample++){
        ResultTree_sampleID[isample] = -2;
        ResultTree_yield[isample] = -2;
        ResultTree_stat_error[isample] = -2;
    }    
    int ResultTree_RunMode = -1;
    ResultTree->SetBranchAddress("ResultTree_sampleSize",&ResultTree_sampleSize);
    ResultTree->SetBranchAddress("ResultTree_plot",&ResultTree_plot);
    ResultTree->SetBranchAddress("ResultTree_sampleID",&ResultTree_sampleID[0]);
    ResultTree->SetBranchAddress("ResultTree_yield",&ResultTree_yield[0]);
    ResultTree->SetBranchAddress("ResultTree_stat_error",&ResultTree_stat_error[0]);
    ResultTree->SetBranchAddress("ResultTree_RunMode",&ResultTree_RunMode);

    float Yields[SignalSamples[1]-SignalSamples[0]+2][3];   // [MC + n signal] [ minus, Normal, plus]
    for(int i=0;i<SignalSamples[1]-SignalSamples[0]+2;i++)
        for(int j=0;j<3;j++) Yields[i][j] = 0.0;

    for(int itree=0;itree<ResultTree->GetEntries();itree++){
        ResultTree->GetEntry(itree);
        if(ResultTree_plot!=RegionIndex) continue;

        if ( ResultTree_RunMode == RunModes_[0]){
            Yields[0][0] = ResultTree_yield[0];
            for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
                Yields[i][0] = ResultTree_yield[SignalSamples[0]+i];
            }
        }else if ( ResultTree_RunMode == Normal){
            Yields[0][1] = ResultTree_yield[0];
            for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
                Yields[i][1] = ResultTree_yield[SignalSamples[0]+i];
            }
        }else if ( ResultTree_RunMode == RunModes_[1]){
            Yields[0][2] = ResultTree_yield[0];
            for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
                Yields[i][2] = ResultTree_yield[SignalSamples[0]+i];
            }
        }
    }

    for(int i=0;i<SignalSamples[1]-SignalSamples[0]+2;i++){
        unc[i] = fabs(Yields[i][0]-Yields[i][1])/Yields[i][1]*100.;
        if(fabs(Yields[i][2]-Yields[i][1])/Yields[i][1]*100. > unc[i]) 
            unc[i] = fabs(Yields[i][2]-Yields[i][1])/Yields[i][1]*100.;
    }
}

void Stat(TChain *ResultTree, int RegionIndex, float unc[]){
    int ResultTree_plot = -1;
    int ResultTree_sampleSize = samples_order_size + 1; // for total MC
    int ResultTree_sampleID[ResultTree_sampleSize];
    double ResultTree_yield[ResultTree_sampleSize];
    double ResultTree_stat_error[ResultTree_sampleSize];
    for(int isample=0;isample<ResultTree_sampleSize;isample++){
        ResultTree_sampleID[isample] = -2;
        ResultTree_yield[isample] = -2;
        ResultTree_stat_error[isample] = -2;
    }    
    int ResultTree_RunMode = -1;
    ResultTree->SetBranchAddress("ResultTree_sampleSize",&ResultTree_sampleSize);
    ResultTree->SetBranchAddress("ResultTree_plot",&ResultTree_plot);
    ResultTree->SetBranchAddress("ResultTree_sampleID",&ResultTree_sampleID[0]);
    ResultTree->SetBranchAddress("ResultTree_yield",&ResultTree_yield[0]);
    ResultTree->SetBranchAddress("ResultTree_stat_error",&ResultTree_stat_error[0]);
    ResultTree->SetBranchAddress("ResultTree_RunMode",&ResultTree_RunMode);

    float Yields[SignalSamples[1]-SignalSamples[0]+2];   // [MC + n signal] 
    float Errors[SignalSamples[1]-SignalSamples[0]+2];   // [MC + n signal] 
    for(int i=0;i<SignalSamples[1]-SignalSamples[0]+2;i++)
    {
        Yields[i] = 0.0;
        Errors[i] = 0.0;
    }

    for(int itree=0;itree<ResultTree->GetEntries();itree++){
        ResultTree->GetEntry(itree);
        if(ResultTree_plot!=RegionIndex) continue;

        if ( ResultTree_RunMode == Normal){
            Yields[0] = ResultTree_yield[0];
            Errors[0] = ResultTree_stat_error[0];
            for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
                Yields[i] = ResultTree_yield[SignalSamples[0]+i];
                Errors[i] = ResultTree_stat_error[SignalSamples[0]+i];
            }
        }
    }

    for(int i=0;i<SignalSamples[1]-SignalSamples[0]+2;i++){
        unc[i] = Errors[i]/Yields[i]*100.;
    }
}

void Xsec(TChain *ResultTree, int RegionIndex, float unc[]){
    int ResultTree_plot = -1;
    int ResultTree_sampleSize = samples_order_size + 1; // for total MC
    int ResultTree_sampleID[ResultTree_sampleSize];
    double ResultTree_yield[ResultTree_sampleSize];
    double ResultTree_stat_error[ResultTree_sampleSize];
    for(int isample=0;isample<ResultTree_sampleSize;isample++){
        ResultTree_sampleID[isample] = -2;
        ResultTree_yield[isample] = -2;
        ResultTree_stat_error[isample] = -2;
    }    
    int ResultTree_RunMode = -1;
    ResultTree->SetBranchAddress("ResultTree_sampleSize",&ResultTree_sampleSize);
    ResultTree->SetBranchAddress("ResultTree_plot",&ResultTree_plot);
    ResultTree->SetBranchAddress("ResultTree_sampleID",&ResultTree_sampleID[0]);
    ResultTree->SetBranchAddress("ResultTree_yield",&ResultTree_yield[0]);
    ResultTree->SetBranchAddress("ResultTree_stat_error",&ResultTree_stat_error[0]);
    ResultTree->SetBranchAddress("ResultTree_RunMode",&ResultTree_RunMode);

    float Yields[SignalSamples[1]-SignalSamples[0]+2][3];   // [MC + n signal] [ minus, Normal, plus]
    for(int i=0;i<SignalSamples[1]-SignalSamples[0]+2;i++)
        for(int j=0;j<3;j++) Yields[i][j] = 0.0;

    // [n MC] [ minus, Normal, plus]
    float SubYields[TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1-
        TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1+1][3];   
    for(int i=0;i<TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1-
            TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1+1;i++)
        for(int j=0;j<3;j++) SubYields[i][j] = 0.0;

    for(int itree=0;itree<ResultTree->GetEntries();itree++){
        ResultTree->GetEntry(itree);
        if(ResultTree_plot!=RegionIndex) continue;

        if ( ResultTree_RunMode == UncXsecMinus){
            Yields[0][0] = ResultTree_yield[0];
            for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
                Yields[i][0] = ResultTree_yield[SignalSamples[0]+i];
            }
            for(int i=0;i<TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1-
                    TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1+1;i++)
                SubYields[i][0] = 
                    ResultTree_yield[TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1
                    +i+1];
        }else if ( ResultTree_RunMode == Normal){
            Yields[0][1] = ResultTree_yield[0];
            for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
                Yields[i][1] = ResultTree_yield[SignalSamples[0]+i];
            }
            for(int i=0;i<TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1-
                    TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1+1;i++)
                SubYields[i][1] = 
                    ResultTree_yield[TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1
                    +i+1];
        }else if ( ResultTree_RunMode == UncXsecPlus){
            Yields[0][2] = ResultTree_yield[0];
            for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
                Yields[i][2] = ResultTree_yield[SignalSamples[0]+i];
            }
            for(int i=0;i<TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1-
                    TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1+1;i++)
                SubYields[i][2] = 
                    ResultTree_yield[TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1
                    +i+1];
        }
    }

    float totalMC_unc = 0.0;
    for(int i=0;i<TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1-
            TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1+1;i++){
        if(fabs(SubYields[i][2]-SubYields[i][1])/Yields[0][1]*100. > 
                fabs(SubYields[i][0]-SubYields[i][1])/Yields[0][1]*100.) {
            totalMC_unc += 
                (fabs(SubYields[i][2]-SubYields[i][1])/Yields[0][1]*100.)*
                (fabs(SubYields[i][2]-SubYields[i][1])/Yields[0][1]*100.);
        }else{
            totalMC_unc += 
                (fabs(SubYields[i][0]-SubYields[i][1])/Yields[0][1]*100.)*
                (fabs(SubYields[i][0]-SubYields[i][1])/Yields[0][1]*100.);
        }
    }
    unc[0] = sqrt(totalMC_unc);
    for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
        unc[i] = fabs(Yields[i][0]-Yields[i][1])/Yields[i][1]*100.;
        if(fabs(Yields[i][2]-Yields[i][1])/Yields[i][1]*100. > unc[i]) 
            unc[i] = fabs(Yields[i][2]-Yields[i][1])/Yields[i][1]*100.;
    }
}

void Scale_Match(TChain *ResultTree, int RegionIndex, float unc[], int UNC_MODE){
    int RunModes_[2] = {0,0};
    int RunModesDY_[2] = {0,0};
    if(UNC_MODE == _Scale){
        RunModes_[0] = TTJets_scaledown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
        RunModes_[1] = TTJets_scaleup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
        RunModesDY_[0] = DYJetsToLL_M_50_scaledown_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
        RunModesDY_[1] = DYJetsToLL_M_50_scaleup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    }else if(UNC_MODE == _Matching){
        RunModes_[0] = TTJets_matchingdown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
        RunModes_[1] = TTJets_matchingup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
        RunModesDY_[0] = DYJetsToLL_M_50_matchingdown_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1;
        RunModesDY_[1] = DYJetsToLL_M_50_matchingup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    }else{
        std::cout<<"[ERROR] This function is not for UNC_MODE of "<<UNC_MODE<<std::endl;
        exit(0);
    }
    int ResultTree_plot = -1;
    int ResultTree_sampleSize = samples_order_size + 1; // for total MC
    int ResultTree_sampleID[ResultTree_sampleSize];
    double ResultTree_yield[ResultTree_sampleSize];
    double ResultTree_stat_error[ResultTree_sampleSize];
    for(int isample=0;isample<ResultTree_sampleSize;isample++){
        ResultTree_sampleID[isample] = -2;
        ResultTree_yield[isample] = -2;
        ResultTree_stat_error[isample] = -2;
    }    
    int ResultTree_RunMode = -1;
    ResultTree->SetBranchAddress("ResultTree_sampleSize",&ResultTree_sampleSize);
    ResultTree->SetBranchAddress("ResultTree_plot",&ResultTree_plot);
    ResultTree->SetBranchAddress("ResultTree_sampleID",&ResultTree_sampleID[0]);
    ResultTree->SetBranchAddress("ResultTree_yield",&ResultTree_yield[0]);
    ResultTree->SetBranchAddress("ResultTree_stat_error",&ResultTree_stat_error[0]);
    ResultTree->SetBranchAddress("ResultTree_RunMode",&ResultTree_RunMode);

    float TotalMC_ = 0.0;
    float Yields[3] = {0.,0.,0.};   // [ minus, Normal, plus]
    float YieldsDY[3] = {0.,0.,0.};   // [ minus, Normal, plus]

    for(int itree=0;itree<ResultTree->GetEntries();itree++){
        ResultTree->GetEntry(itree);
        if(ResultTree_plot!=RegionIndex) continue;

        if(ResultTree_RunMode == UncQsquare ){
            Yields[0] = ResultTree_yield[(int)(RunModes_[0]+1)];
            Yields[1] = 
                ResultTree_yield[TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1+1];
            Yields[2] = ResultTree_yield[(int)(RunModes_[1]+1)];

            YieldsDY[0] = ResultTree_yield[(int)(RunModesDY_[0]+1)];
            YieldsDY[1] = 
                ResultTree_yield[DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1+1];
            YieldsDY[2] = ResultTree_yield[(int)(RunModesDY_[1]+1)];
            TotalMC_ = ResultTree_yield[0];
            //std::cout<<" ( "<<Yields[0]<<"-"<<RunModes_[0] <<"- , "<<Yields[1]<< " , "<<Yields[2]
            //    <<"-"<< RunModes_[1]<<"- , "<<TotalMC_<<std::endl;
            //std::cout<<" ( "<<YieldsDY[0]<<"-"<<RunModesDY_[0] <<"- , "<<YieldsDY[1]<< " , "<<YieldsDY[2]
            //    <<"-"<< RunModesDY_[1]<<"- , "<<TotalMC_<<std::endl;
        }
    }

    unc[0] = fabs(Yields[0]-Yields[1])/TotalMC_*100.;
    if(fabs(Yields[2]-Yields[1])/TotalMC_*100.>unc[0])
        unc[0] = fabs(Yields[2]-Yields[1])/TotalMC_*100.;

    float unc_ = fabs(YieldsDY[0]-YieldsDY[1])/TotalMC_*100.;
    if(fabs(YieldsDY[2]-YieldsDY[1])/TotalMC_*100.>unc_)
        unc_ = fabs(YieldsDY[2]-YieldsDY[1])/TotalMC_*100.;

    unc[0] = sqrt(unc[0]*unc[0] + unc_*unc_);
}

void Lep_Pho_ID(TChain *ResultTree, int RegionIndexB4, int RegionIndexaf, float unc[]){
    int ResultTree_plot = -1;
    int ResultTree_sampleSize = samples_order_size + 1; // for total MC
    int ResultTree_sampleID[ResultTree_sampleSize];
    double ResultTree_yield[ResultTree_sampleSize];
    double ResultTree_stat_error[ResultTree_sampleSize];
    for(int isample=0;isample<ResultTree_sampleSize;isample++){
        ResultTree_sampleID[isample] = -2;
        ResultTree_yield[isample] = -2;
        ResultTree_stat_error[isample] = -2;
    }    
    int ResultTree_RunMode = -1;
    ResultTree->SetBranchAddress("ResultTree_sampleSize",&ResultTree_sampleSize);
    ResultTree->SetBranchAddress("ResultTree_plot",&ResultTree_plot);
    ResultTree->SetBranchAddress("ResultTree_sampleID",&ResultTree_sampleID[0]);
    ResultTree->SetBranchAddress("ResultTree_yield",&ResultTree_yield[0]);
    ResultTree->SetBranchAddress("ResultTree_stat_error",&ResultTree_stat_error[0]);
    ResultTree->SetBranchAddress("ResultTree_RunMode",&ResultTree_RunMode);

    float Yields[SignalSamples[1]-SignalSamples[0]+2][3];   // [MC + n signal] [ minus, Normal, plus]
    for(int i=0;i<SignalSamples[1]-SignalSamples[0]+2;i++)
        for(int j=0;j<3;j++) Yields[i][j] = 0.0;

    for(int itree=0;itree<ResultTree->GetEntries();itree++){
        ResultTree->GetEntry(itree);
        if(ResultTree_RunMode!=Normal) continue;

        if ( ResultTree_plot == RegionIndexB4){
            Yields[0][0] = ResultTree_yield[0];
            for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
                Yields[i][0] = ResultTree_yield[SignalSamples[0]+i];
            }
        }else if ( ResultTree_plot == RegionIndexaf){
            Yields[0][1] = ResultTree_yield[0];
            for(int i=1;i<SignalSamples[1]-SignalSamples[0]+2;i++){
                Yields[i][1] = ResultTree_yield[SignalSamples[0]+i];
            }
        }
    }

    for(int i=0;i<SignalSamples[1]-SignalSamples[0]+2;i++){
        unc[i] = fabs(Yields[i][0]-Yields[i][1])/Yields[i][1]*100.;
    }
}

void TrimChar(char *OrignalChar_, char *TrimChar_){
    int Nchar = strlen(OrignalChar_) -6 -3;  // remove Tgamma and GeV
    for(int i=0;i<Nchar;i++)
        TrimChar_[i] = OrignalChar_[6+i];
}

