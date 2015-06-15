#include "TTree.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TFile.h"
#include "interface/format_mini.h"
//#include "interface/format.h"
#include "interface/MCTruth.h"
#include "interface/TriggerBooking.h"
#include "interface/MCTruthChannel.h"
//#include "interface/HitFitInfoBranches.h"
#include "interface/DPHI.h"
#include "interface/DR.h"
#include "interface/JER.h"
//#include "interface/samples_tr_mu.h"
#include "interface/samples_tr_1lepton4jets.h"
#include "interface/ScalePhotonID.h"
#include "interface/ScaleLeptonID.h"
#include "interface/ObjectMCTruthForRatio.h"
#include "interface/IsTTbarDuplicateRemoval.h"
//#include "interface/samples_tr.h"

int RunStatus_ = Normal;   
/*
enum RunStatus{
    Normal,
    UncXsecPlus,
    UncXsecMinus,
    UncPUPlus,
    UncPUMinus,
    UncJESPlus,
    UncJESMinus,
    UncJERPlus,
    UncJETMinus,
    UncQsquare,
    RunStatusSize
};
*/
#include "interface/DRAWSTACK.h"
#include "interface/HistMerge.h"
#include "interface/checkEvt.h"
#include "interface/RecoLooseLeptonSelection.h"
#include "interface/RecoLeptonSelection.h"
#include "interface/RecoLeptonForCleaningSelection.h"
#include "interface/RecoPhotonSelection.h"
#include "interface/RecoPhotonForJetCleaningSelection.h"
#include "interface/RecoLoosePhotonSelection.h"
#include "interface/RecoPhotonSelectionCheckingDR.h"
#include "interface/RecoJetSelection.h"
#include "interface/SolutionOfWNeutrino.h"
#include "interface/ReduceTree.h"
#include "interface/PUReweighting.h"
#include "interface/TopTriggerEfficiencyProvider.cc"
//#include "interface/TMVAClassification_CutsD.class.C"
#include "TChain.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "interface/RecoZbosonFordiPhoton.h"
#include "interface/RecoTstar.h"
#include "interface/RecoFakeTstar.h"
#include "interface/dMRecoTstar.h"
#include <math.h>
#include <string>
#include <iostream>
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "stdlib.h"
#include <fstream>
#include <map>
#include <vector>
#include <time.h>
#include <math.h>
#include "interface/RecoZboson.h"
#include "interface/RecoMoflljj.h"
#include "interface/RecoMZb.h"
#include "interface/PhotonMCTruth.h"
#include "interface/RecoNoCutPhotonSelection.h"
#include "interface/CorrectionOnTstar.h"
#include "interface/MainSysUnc.h"



bool RunForCR_Shape = false;
bool RunForMoreStat = false;

double Weight3D_[50][50][50];
const int Nrebins = 5;
const int NRuns =5;
string RunsName[NRuns] = {"A","B","C","D","ALL"};
TH1D *PileUpReweight[NRuns]; // RunA/B/C/D/all
float LUMINOSITYForRuns[NRuns] = {876.2,4411.7,7055.2,7369.0,19712.1};
float JetPtOrderCuts[NRuns][4] = {
    {55,45,35,20},
    {55,45,35,20},
    {55,45,35,20},
    {55,45,35,20},
    {55,45,35,20}
};
int RunsRanges[NRuns][2] = {
    {190645,193621},
    {193834,196531},
    {198049,203002},
    {203709,208686},
    {190645,208686}
};
TopTriggerEfficiencyProvider *weight_provider[NRuns] ;

TH2D *hRatioData;
TH2D *hRatioMCBG;

const int Nratio = 8;
TH2D *hRatioMCBGpTandEta[Nratio];    // [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
TH2D *hMCBGpTandEtaNumerator[Nratio];    // [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
TH2D *hMCBGpTandEtaDenominator[Nratio];    // [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
TH1D *hMCBGpTNumerator[Nratio];    // [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
TH1D *hMCBGpTDenominator[Nratio];    // [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
TH1D *hMCBGEtaNumerator[Nratio];    // [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
TH1D *hMCBGEtaDenominator[Nratio];    // [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
string RatioNames[Nratio] = {
    "photon",
    "ISR_FSR",
    "gluon",
    "light_quark",
    "heavy_quark",
    "electron",
    "muon",
    "tau"
};

enum plots_order{
    _N1lep4jet,
    _NOfVertexb4PU1l4j,
    _NOfVertexafPU1l4j,
    _NOfVertexafPUTrgeff1l4j,
    _N1p1lep4jet,
    _NOfVertexb4PU1p1l4j,
    _NOfVertexb4PU1p1l4jWithTstarCorrection,
    _NOfVertexafPU1p1l4j,
    _NOfVertexafPUTrgeff1p1l4j,
    _NOfVertexafPUTrgeff1p1l4jafLS,
    _NOfVertexafPUTrgeff1p1l4jafLSPS,
    _NOfVertexafPUTrgeff1p1l4jafLSPSTop,
    _Nj,
    _Np,
    _Ne,
    _Nu,
    _PhotonPt,
    _JetPt,
    _JetPtRunA,
    _JetPtRunB,
    _JetPtRunC,
    _JetPtRunD,
    _Jet1stPt,
    _Jet2ndPt,
    _Jet3thPt,
    _Jet4thPt,
    _ElectronPt,
    _MuonPt,
    _MET,
    _ST, 
    _N1p1lep4jetRunA,
    _N1p1lep4jetRunB,
    _N1p1lep4jetRunC,
    _N1p1lep4jetRunD,
    _N1p1lep5jet,
    _N1lep6jet,
    _N1p2lep4jet,
    _N3lep4jet,
    _N2lep5jet,
    _N2lep5jetRunA,
    _N2lep5jetRunB,
    _N2lep5jetRunC,
    _N2lep5jetRunD,
    _N2p1lep4jet,
    _MtstarChi2_2r1l4jnoLSPS,
    _MtstarChi2_2r1l4jnoPS,
    _MtstarChi2_2r1l4j,
    _MtstarChi2_2r1l4jReweight1,
    _MtstarChi2_2r1l4jReweight2,
    _MtstarChi2_1r1l5j,
    _MtstarChi2_0r1l6j,
    _MtstarChi2_1r2l4j,
    _MtstarChi2_0r3l4j,
    _MtstarChi2_0r2l5j,
    _NjInSR,
    _NpInSR,
    _NeInSR,
    _NuInSR,
    _PhotonPtInSR,
    _JetPtInSR,
    _JetPtRunAInSR,
    _JetPtRunBInSR,
    _JetPtRunCInSR,
    _JetPtRunDInSR,
    _Jet1stPtInSR,
    _Jet2ndPtInSR,
    _Jet3thPtInSR,
    _Jet4thPtInSR,
    _ElectronPtInSR,
    _MuonPtInSR,
    _METInSR,
    _STInSR, 
    _plots_order_size
};

struct plots_struct{
    char PlotName[128];
    char YaxisTitle_[128];
    char XaxisTitle_[128];
    char XUnit_[128];
    int Ndivide_;
    double XminBond_;
    double XmaxBond_;
    int rebin_;
};

struct plots_struct PLOTS[_plots_order_size] = {
    {"_N1lep4jet","Yields","_N1lep4jet","",1,0,1,1},
    {"_NOfVertexb4PU1l4j","Yields","N(Vertices) b4PU : 1l4j","-",60,0,60,1},
    {"_NOfVertexafPU1l4j","Yields","N(Vertices) afPU : 1l4j","-",60,0,60,1},
    {"_NOfVertexafPUTrgeff1l4j","Yields","N(Vertices) afPU+Trg : 1l4j","-",60,0,60,1},
    {"_N1p1lep4jet","Yields","_N1p1lep4jet","",1,0,1,1},
    {"_NOfVertexb4PU1p1l4j","Yields","N(Vertices) b4PU : 1p1l4j","-",60,0,60,1},
    {"_NOfVertexb4PU1p1l4jWithTstarCorrection","Yields","N(Vertices) b4PU : 1p1l4j","-",60,0,60,1},
    {"_NOfVertexafPU1p1l4j","Yields","N(Vertices) afPU : 1p1l4j","-",60,0,60,1},
    {"_NOfVertexafPUTrgeff1p1l4j","Yields","N(Vertices) afPU+Trg : 1p1l4j","-",60,0,60,1},
    {"_NOfVertexafPUTrgeff1p1l4jafLS","Yields","N(Vertices) afPU+Trg+LS : 1p1l4j","-",60,0,60,1},
    {"_NOfVertexafPUTrgeff1p1l4jafLSPS","Yields","N(Vertices) afPU+Trg+LS+PS : 1p1l4j","-",60,0,60,1},
    {"_NOfVertexafPUTrgeff1p1l4jafLSPSTop","Yields","N(Vertices) afPU+Trg+LS+PS+Top : 1p1l4j","-",60,0,60,1},

    {"_Nj","Yields","N(jets)","-",11,0,11,1},
    {"_Np","Yields","N(photons)","-",11,0,11,1},
    {"_Ne","Yields","N(electrons)","-",11,0,11,1},
    {"_Nu","Yields","N(muons)","-",11,0,11,1},
    {"_PhotonPt","Yields","Photon Pt","GeV/c",200,0,200,Nrebins},
    {"_JetPt","Yields","Jet Pt","GeV/c",200,0,200,Nrebins},
    {"_JetPtRunA","Yields","Jet Pt(RunA)","GeV/c",200,0,200,Nrebins},
    {"_JetPtRunB","Yields","Jet Pt(RunB)","GeV/c",200,0,200,Nrebins},
    {"_JetPtRunC","Yields","Jet Pt(RunC)","GeV/c",200,0,200,Nrebins},
    {"_JetPtRunD","Yields","Jet Pt(RunD)","GeV/c",200,0,200,Nrebins},
    {"_Jet1stPt","Yields","Jet 1st Pt","GeV/c",200,0,200,Nrebins},
    {"_Jet2ndPt","Yields","Jet 2nd Pt","GeV/c",200,0,200,Nrebins},
    {"_Jet3thPt","Yields","Jet 3th Pt","GeV/c",200,0,200,Nrebins},
    {"_Jet4thPt","Yields","Jet 4th Pt","GeV/c",200,0,200,Nrebins},
    {"_ElectronPt","Yields","Electron Pt","GeV/c",200,0,200,Nrebins},
    {"_MuonPt","Yields","Muon Pt","GeV/c",200,0,200,Nrebins},
    {"_MET","Yields","MET","GeV/c",200,0,200,Nrebins},
    {"_ST","Yields","ST","GeV/c",1800,200,2000,50},

    {"_N1p1lep4jetRunA","Yields","_N1p1lep4jetRunA","",1,0,1,1},
    {"_N1p1lep4jetRunB","Yields","_N1p1lep4jetRunB","",1,0,1,1},
    {"_N1p1lep4jetRunC","Yields","_N1p1lep4jetRunC","",1,0,1,1},
    {"_N1p1lep4jetRunD","Yields","_N1p1lep4jetRunD","",1,0,1,1},
    {"_N1p1lep5jet","Yields","_N1p1lep5jet","",1,0,1,1},
    {"_N1lep6jet","Yields","_N1lep6jet","",1,0,1,1},
    {"_N1p2lep4jet","Yields","_N1p2lep4jet","",1,0,1,1},
    {"_N3lep4jet","Yields","_N3lep4jet","",1,0,1,1},
    {"_N2lep5jet","Yields","_N2lep5jet","",1,0,1,1},
    {"_N2lep5jetRunA","Yields","_N2lep5jetRunA","",1,0,1,1},
    {"_N2lep5jetRunB","Yields","_N2lep5jetRunB","",1,0,1,1},
    {"_N2lep5jetRunC","Yields","_N2lep5jetRunC","",1,0,1,1},
    {"_N2lep5jetRunD","Yields","_N2lep5jetRunD","",1,0,1,1},
    {"_N2p1lep4jet","Yields","_N2p1lep4jet","",1,0,1,1},
    {"_MtstarChi2_2r1l4jnoLSPS","Yields","M(t*:2r1l4j) w/o LS+PS","GeV/c^{2}",1500,50,1550,150},
    {"_MtstarChi2_2r1l4jnoPS","Yields","M(t*:2r1l4j) w/o PS","GeV/c^{2}",1500,50,1550,150},
    {"_MtstarChi2_2r1l4j","Yields","M(t*:2r1l4j)","GeV/c^{2}",1500,50,1550,150},
    {"_MtstarChi2_2r1l4jReweight1","Yields","M(t*:2r1l4j,reweight1)","GeV/c^{2}",1500,50,1550,150},
    {"_MtstarChi2_2r1l4jReweight2","Yields","M(t*:2r1l4j,reweight2)","GeV/c^{2}",1500,50,1550,150},
    {"_MtstarChi2_1r1l5j","Yields","M(t*:1r1l5j)","GeV/c^{2}",1500,50,1550,150},
    {"_MtstarChi2_0r1l6j","Yields","M(t*:0r1l6j)","GeV/c^{2}",1500,50,1550,150},
    {"_MtstarChi2_1r2l4j","Yields","M(t*:1r2l4j)","GeV/c^{2}",1500,50,1550,150},
    {"_MtstarChi2_0r3l4j","Yields","M(t*:0r3l4j)","GeV/c^{2}",1500,50,1550,150},
    {"_MtstarChi2_0r2l5j","Yields","M(t*:0r2l5j)","GeV/c^{2}",1500,50,1550,150},

    {"_NjInSR","Yields","N(jets)","-",11,0,11,1},
    {"_NpInSR","Yields","N(photons)","-",11,0,11,1},
    {"_NeInSR","Yields","N(electrons)","-",11,0,11,1},
    {"_NuInSR","Yields","N(muons)","-",11,0,11,1},
    {"_PhotonPtInSR","Yields","Photon Pt","GeV/c",200,0,200,20},
    {"_JetPtInSR","Yields","Jet Pt","GeV/c",200,0,200,20},
    {"_JetPtRunAInSR","Yields","Jet Pt(RunA)","GeV/c",200,0,200,20},
    {"_JetPtRunBInSR","Yields","Jet Pt(RunB)","GeV/c",200,0,200,20},
    {"_JetPtRunCInSR","Yields","Jet Pt(RunC)","GeV/c",200,0,200,20},
    {"_JetPtRunDInSR","Yields","Jet Pt(RunD)","GeV/c",200,0,200,20},
    {"_Jet1stPtInSR","Yields","Jet 1st Pt","GeV/c",200,0,200,20},
    {"_Jet2ndPtInSR","Yields","Jet 2nd Pt","GeV/c",200,0,200,20},
    {"_Jet3thPtInSR","Yields","Jet 3th Pt","GeV/c",200,0,200,20},
    {"_Jet4thPtInSR","Yields","Jet 4th Pt","GeV/c",200,0,200,20},
    {"_ElectronPtInSR","Yields","Electron Pt","GeV/c",200,0,200,20},
    {"_MuonPtInSR","Yields","Muon Pt","GeV/c",200,0,200,20},
    {"_METInSR","Yields","MET","GeV/c",200,0,200,20},
    {"_STInSR","Yields","ST","GeV/c",1800,200,2000,150}
};

enum plots2D_order{
    _NjetNlepMC,
    _NjetNlepData,
    _NjetNlepRunA,
    _NjetNlepRunB,
    _NjetNlepRunC,
    _NjetNlepRunD,
    _1tp1lep4jet,
    _1lp1lep4jet,
    _1tp1lep4jetExcludingB1,
    _1lp1lep4jetExcludingB1,
    _ABCDdMST,
    _plots2D_order_size
};

struct plots2D_struct{
    char PlotName[128];
    char YaxisTitle_[128];
    char XaxisTitle_[128];
    char YUnit_[128];
    int YNdivide_;
    double YminBond_;
    double YmaxBond_;
    int Yrebin_;
    char XUnit_[128];
    int XNdivide_;
    double XminBond_;
    double XmaxBond_;
    int Xrebin_;
};
struct plots2D_struct PLOTS2D[_plots2D_order_size] = { 
    {"_NjetNlepMC","Njet(MC)","Nlep(MC)","-",10,0,10,1,"-",10,0,10,1},
    {"_NjetNlepData","Njet(Data)","Nlep(Data)","-",10,0,10,1,"-",10,0,10,1},
    {"_NjetNlepRunA","Njet(RunA)","Nlep(RunA)","-",10,0,10,1,"-",10,0,10,1},
    {"_NjetNlepRunB","Njet(RunB)","Nlep(RunB)","-",10,0,10,1,"-",10,0,10,1},
    {"_NjetNlepRunC","Njet(RunC)","Nlep(RunC)","-",10,0,10,1,"-",10,0,10,1},
    {"_NjetNlepRunD","Njet(RunD)","Nlep(RunD)","-",10,0,10,1,"-",10,0,10,1},
    {"_1tp1lep4jet","pT(#gamma)","#eta(#gamma)","-",200,0,200,1,"-",100,-2.5,2.5,2},
    {"_1lp1lep4jet","pT(#gamma)","#eta(#gamma)","-",200,0,200,1,"-",100,-2.5,2.5,2},
    {"_1tp1lep4jetExcludingB1","pT(#gamma)","#eta(#gamma)","-",200,0,200,1,"-",100,-2.5,2.5,2},
    {"_1lp1lep4jetExcludingB1","pT(#gamma)","#eta(#gamma)","-",200,0,200,1,"-",100,-2.5,2.5,2},
    {"_ABCDdMST","abs(M1-M2)/(M1+M2)","ST","-",20,0,1,1,"GeV",3000,0,3000,100}
};

TH1F *hMainSysUnc_[_plots_order_size];
/*
   For Excited Quark : 
   (lepb_label+Wlep+gluon1_label) = (hadb_label+hadw1_label+hadw2_label+gluon2_label)
   (11 + Wlep + 16 ) = (12 + 13 + 14 + 17)
 */

void analyze(int sample_idx,TH1F *ExcitedQuarkPlots[], TH2D *ExcitedQuark2DPlots[], double weight_, int MCTruthChannels[])
{


    TChain *root = new TChain("bprimeKit/root");
    root->Add(SAMPLE[sample_idx].filename);

    EvtInfoBranches EvtInfo;
    EvtInfo.Register(root);
    GenInfoBranches GenInfo;
    GenInfo.Register(root);
    LepInfoBranches LepInfo;
    LepInfo.Register(root,"PFLepInfo");
    JetInfoBranches JetInfo;
    JetInfo.Register(root,"PFJetInfo");
    VertexInfoBranches VertexInfo;
    VertexInfo.Register(root);
    PhotonInfoBranches PhotonInfo;
    PhotonInfo.Register(root);

    ReduceTree(root);

    double weight = weight_;

    char buffer[128];
    sprintf(buffer,"tree_%s%s",SAMPLE[sample_idx].tag,RunStatusNames[RunStatus_].c_str());
    TTree *MCTemplatesTree = new TTree(buffer,buffer);

    int sample_ID=sample_idx;
    int Nr_template=-1;
    int Nj_template=-1;
    int Nl_template=-1;
    int Nlloose_template=-1;
    int Nrloose_template=-1;
    float weight_template=-1;
    float Mass_template=-1;
    float category_template=-1;
    int MCTruth1 = -1;
    int MCTruth2 = -1;
    int MCTruth3 = -1;
    int MCTruth_tmp[3] = {-1,-1,-1};

    MCTemplatesTree->Branch("sampleID",&sample_ID,"sampleID/I");
    MCTemplatesTree->Branch("weight",&weight_template,"weight/F");
    MCTemplatesTree->Branch("Mass",&Mass_template,"Mass/F");
    MCTemplatesTree->Branch("category",&category_template,"category/F");
    MCTemplatesTree->Branch("Nr",&Nr_template,"Nr/I");
    MCTemplatesTree->Branch("Nj",&Nj_template,"Nj/I");
    MCTemplatesTree->Branch("Nl",&Nl_template,"Nl/I");
    MCTemplatesTree->Branch("Nl_loose",&Nlloose_template,"Nl_loose/I");
    MCTemplatesTree->Branch("Nr_loose",&Nrloose_template,"Nr_loose/I");
    MCTemplatesTree->Branch("MCTruth1",&MCTruth1,"MCTruth1/I");
    MCTemplatesTree->Branch("MCTruth2",&MCTruth2,"MCTruth2/I");
    MCTemplatesTree->Branch("MCTruth3",&MCTruth3,"MCTruth3/I");

    int RunDependies[NRuns] = {
        _N1p1lep4jetRunA,
        _N1p1lep4jetRunB,
        _N1p1lep4jetRunC,
        _N1p1lep4jetRunD,
        -1  // no use for all run mode
    };
    int JetPtRunIdx[NRuns] = {
        _JetPtRunA,
        _JetPtRunB,
        _JetPtRunC,
        _JetPtRunD,
        -1  // no use for all run mode
    };



    for(int i=0;i<_plots_order_size;i++) ExcitedQuarkPlots[i]->Sumw2();

    map< pair<int, int>, int > evtlist[NRuns];
    int nevents_total = root->GetEntries();	

    int PhoMode_ = 1;
    int LepMode_ = 1;
    int TopPtMode_ = 1;

    switch(RunStatus_){
        case UncPhoIDPlus:
            PhoMode_ = 2;
            break;
        case UncPhoIDMinus:
            PhoMode_ = 0;
            break;
        case UncLepIDPlus:
            LepMode_ = 2;
            break;
        case UncLepIDMinus:
            LepMode_ = 0;
            break;
        case UncTopPtPlus:
            TopPtMode_ = 2;
            break;
        case UncTopPtMinus:
            TopPtMode_ = 0;
            break;
        default:
            PhoMode_ = 1;
            LepMode_ = 1;
            TopPtMode_ = 1;
            break;
    }

    // PDF information
    FILE *pdffileCR;
    FILE *pdffileSR;

    sprintf(buffer,"PDF/%s_%s.txt",SAMPLE[sample_idx].tag,PLOTS[_ST].PlotName);
    pdffileCR = fopen(buffer,"w");
    sprintf(buffer,"PDF/%s_%s.txt",SAMPLE[sample_idx].tag,PLOTS[_STInSR].PlotName);
    pdffileSR = fopen(buffer,"w");


    for(int entry=0;entry<root->GetEntries();entry++) {
        if ((entry%10000) == 0)	printf("Loading event #%d of %d.\n",entry,nevents_total);
        root->GetEntry(entry);

        // remove ttbar events overlapping with ttG sample
        if((sample_idx >= TTJets_matchingdown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1 &&
                    sample_idx<=TTJets_scaleup_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1) || 
                sample_idx == TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1)
            if(IsTTbarDuplicateRemoval( GenInfo)) continue;

        char DATA_TAG[128] ="data";
        if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
            // Smearing MC Jet's pT (https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#2012)
            int JERMODE = 0;
            for(int nj=0;nj<JetInfo.Size;nj++){
                if(JetInfo.GenJetPt[nj]>1){
                    switch(RunStatus_){
                        case UncJERPlus:
                            JERMODE = 1;
                            break;
                        case UncJERMinus:
                            JERMODE = -1;
                            break;
                        default:
                            JERMODE = 0;
                            break;
                    }
                    JetInfo.Pt[nj] = 
                        max(0.,(double)(JetInfo.GenJetPt[nj] 
                                    + JER(JetInfo.Eta[nj],JERMODE)*(JetInfo.Pt[nj]-JetInfo.GenJetPt[nj])));
                }
            }
            // For JES
            int JESMODE = 0;
            for(int nj=0;nj<JetInfo.Size;nj++){
                switch(RunStatus_){
                    case UncJESPlus:
                        JESMODE = 1;
                        break;
                    case UncJESMinus:
                        JESMODE = -1;
                        break;
                    default:
                        JESMODE = 0;
                        break;
                }
                JetInfo.Pt[nj] = JetInfo.Pt[nj]*(1+JESMODE*fabs(JetInfo.Unc[nj]));
            }

        }

        // some plots for RunA/B/C/D separately
        for(int irun_=0;irun_<NRuns-1;irun_++){

            if(!strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
                if(!(EvtInfo.RunNo>=RunsRanges[irun_][0]&&EvtInfo.RunNo<=RunsRanges[irun_][1])) continue;
                if(!isGoodEvt(EvtInfo.RunNo,EvtInfo.LumiNo)) continue;

                // remove duplicate event
                map< pair<int, int> , int>::iterator evtitr;
                evtitr = evtlist[irun_].find( pair<int, int>(EvtInfo.RunNo, EvtInfo.EvtNo) );
                if( evtitr == evtlist[irun_].end() )
                    evtlist[irun_].insert( pair<pair<int, int>, int>(pair<int, int>(EvtInfo.RunNo, EvtInfo.EvtNo), 1));
                else
                    continue;
            }

            // Object selection -- start --
            int NMuons = 0;
            int NElectrons = 0;
            int NLeptons = 0;
            int NMuonsC = 0;
            int NElectronsC = 0;
            int NLeptonsC = 0;
            int NMuonsV = 0;
            int NElectronsV = 0;
            int NLeptonsV = 0;
            int NPhotons = 0;
            int NPhotonsNoISO = 0;
            int NJets = 0;
            int M_Index[MAX_LEPTONS/2];
            int E_Index[MAX_LEPTONS/2];
            int L_Index[MAX_LEPTONS/2*2];
            int MC_Index[MAX_LEPTONS/2];
            int EC_Index[MAX_LEPTONS/2];
            int LC_Index[MAX_LEPTONS/2*2];
            int MV_Index[MAX_LEPTONS/2];
            int EV_Index[MAX_LEPTONS/2];
            int LV_Index[MAX_LEPTONS/2*2];
            int P_Index[MAX_PHOTONS];
            int PnoISO_Index[MAX_PHOTONS];
            int J_Index[MAX_JETS];
            for(int m=0;m<MAX_LEPTONS/2;m++) M_Index[m] = 0;
            for(int e=0;e<MAX_LEPTONS/2;e++) E_Index[e] = 0;
            for(int l=0;l<MAX_LEPTONS/2*2;l++) L_Index[l] = 0;

            for(int m=0;m<MAX_LEPTONS/2;m++) MC_Index[m] = 0;
            for(int e=0;e<MAX_LEPTONS/2;e++) EC_Index[e] = 0;
            for(int l=0;l<MAX_LEPTONS/2*2;l++) LC_Index[l] = 0;
            for(int m=0;m<MAX_LEPTONS/2;m++) MV_Index[m] = 0;
            for(int e=0;e<MAX_LEPTONS/2;e++) EV_Index[e] = 0;
            for(int l=0;l<MAX_LEPTONS/2*2;l++) LV_Index[l] = 0;

            for(int p=0;p<MAX_PHOTONS;p++) P_Index[p] = 0;
            for(int p=0;p<MAX_PHOTONS;p++) PnoISO_Index[p] = 0;
            for(int j=0;j<MAX_JETS;j++) J_Index[j] = 0;
            RecoLooseLeptonSelection(EvtInfo,LepInfo,NMuonsV,MV_Index,NElectronsV,EV_Index,NLeptonsV,LV_Index);
            RecoLeptonSelection(EvtInfo,LepInfo,NMuons,M_Index,NElectrons,E_Index,NLeptons,L_Index);
            RecoLeptonForCleaningSelection(EvtInfo,LepInfo,NMuonsC,MC_Index,NElectronsC,EC_Index,NLeptonsC,LC_Index);

            UsingPhoISO = 10000000000.0;    // meaning no ISO reqirement
            //RecoPhotonSelection(LepInfo,PhotonInfo,NMuonsV,MV_Index,NElectronsV,EV_Index,NPhotonsNoISO,PnoISO_Index,EvtInfo.RhoPU[0],P_LOOSE);
            RecoPhotonSelection(LepInfo,PhotonInfo,NMuonsC,MC_Index,NElectronsC,EC_Index,NPhotonsNoISO,PnoISO_Index,EvtInfo.RhoPU[0],P_LOOSE);

            UsingPhoISO = 1.0;
            //RecoPhotonSelection(LepInfo,PhotonInfo,NMuonsV,MV_Index,NElectronsV,EV_Index,NPhotons,P_Index,EvtInfo.RhoPU[0],P_LOOSE);
            RecoPhotonSelection(LepInfo,PhotonInfo,NMuonsC,MC_Index,NElectronsC,EC_Index,NPhotons,P_Index,EvtInfo.RhoPU[0],P_LOOSE);
            // there is a difference using standard lepton (muon) for photon cleaning
            //RecoPhotonSelection(LepInfo,PhotonInfo,NMuons,M_Index,NElectrons,E_Index,NPhotons,P_Index,EvtInfo.RhoPU[0],P_LOOSE);

            RecoJetSelection(LepInfo,JetInfo, PhotonInfo,NMuonsC,MC_Index,NElectronsC,EC_Index,NJets,J_Index,NPhotons,P_Index);

            // number of vertex -- start --
            int NOfVertex_=0;
            int V_Index[MAX_Vertices];
            for(int v=0;v<MAX_Vertices;v++) V_Index[v] = 0;
            for(int j=0;j<VertexInfo.Size;j++) {
                if(VertexInfo.isValid[j]!=1) continue; 
                if(VertexInfo.Type[j]!=0) continue; 
                if(VertexInfo.Ndof[j]<=4) continue; 
                if(fabs(VertexInfo.Rho[j])>=2) continue; 
                if(fabs(VertexInfo.z[j])>=24) continue; 
                V_Index[NOfVertex_] = j;
                NOfVertex_++;
            }
            if(NOfVertex_<1) continue;

            if((NLeptons>=1&&NJets>=4)) {
                if(JetInfo.Pt[J_Index[0]]>JetPtOrderCuts[irun_][0]) 
                if(JetInfo.Pt[J_Index[1]]>JetPtOrderCuts[irun_][1]) 
                    if(JetInfo.Pt[J_Index[2]]>JetPtOrderCuts[irun_][2]){ 
                        bool IsMuon =0;
                        if (LepInfo.LeptonType[L_Index[0]]==13) IsMuon =1;
                        if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG)) {

                            float WeightCorrectionOnTstar = 1.0;
                            if(sample_idx>=Tgamma500GeV&&sample_idx<=TGluonTgamma1500GeV){
                                int MCTruthChannel_ = MCTruthChannel(GenInfo);
                                WeightCorrectionOnTstar = CorrectionOnTstar(sample_idx,MCTruthChannel_);
                            }

                            if((EvtInfo.TrueIT[0])>=60){
                                std::cout<<"[PU WARNING] N(vertex = "<<EvtInfo.TrueIT[0]<<") >= 60"<<std::endl;
                            }else{
                                weight = WeightCorrectionOnTstar*weight_*PileUpReweight[irun_]->GetBinContent(PileUpReweight[irun_]->GetXaxis()->FindBin(EvtInfo.TrueIT[0]))*LUMINOSITYForRuns[irun_]/LUMINOSITY;
                            }

                            // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Trigger_Efficiencies
                            std::vector<double> WeightTopTriggerEfficiency = weight_provider[irun_]->get_weight(
                                    LepInfo.Pt[L_Index[0]],LepInfo.Eta[L_Index[0]],
                                    JetInfo.Pt[J_Index[3]],JetInfo.Eta[J_Index[3]],
                                    NOfVertex_,
                                    NJets,
                                    IsMuon,
                                    TopTriggerEfficiencyProvider::NOMINAL);
                            double uncOnTrig = 0.;
                            if(RunStatus_ == UncTrigPlus) {
                                uncOnTrig = 1.0;
                                //uncOnTrig = 0.01;
                            }else if (RunStatus_ == UncTrigMinus){
                                uncOnTrig = -1.0;
                                //uncOnTrig = -0.01;
                            }
                            weight *= (WeightTopTriggerEfficiency[0] + uncOnTrig*WeightTopTriggerEfficiency[1]);
                            //weight *= (WeightTopTriggerEfficiency[0] + uncOnTrig);
                            if(NPhotons>=1&&NLeptonsV==1){
                                float scaleFactor_ =  
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]],LepMode_)*
                                    ScalePhotonID(PhotonInfo.Pt[P_Index[0]], PhotonInfo.Eta[P_Index[0]], P_LOOSE, PhoMode_);

                                float TopPtReweighting = 1.0;
                                if (sample_idx == 
                                TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1){
                                    int idxes[2] = {-1,-1};
                                    for(int iTop=0;iTop<GenInfo.Size;iTop++){
                                        if(GenInfo.PdgID[iTop] == 6) idxes[0] = iTop; 
                                        if(GenInfo.PdgID[iTop] == -6) idxes[1] = iTop; 
                                    }   
                                    if(idxes[0]!=-1&&idxes[1]!=-1){
                                        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
                                        TopPtReweighting = 
                                        sqrt(exp(0.159-0.00141*GenInfo.Pt[idxes[0]])*exp(0.159-0.00141*GenInfo.Pt[idxes[1]]));

                                        TopPtReweighting = pow(TopPtReweighting,TopPtMode_);
                                    }else{
                                        std::cout<<"[WARNING] no Gen top"<<std::endl;
                                    }

                                }
                                weight = weight*TopPtReweighting;
                            }
                        }

                        ExcitedQuark2DPlots[_NjetNlepRunA+irun_]->Fill(NJets,NLeptons,weight);

                        // loose regions
                        if(NPhotons>=1&&NLeptonsV==1){
                            float scaleFactor_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG))
                                scaleFactor_= 
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_)*
                                    ScalePhotonID(PhotonInfo.Pt[P_Index[0]], PhotonInfo.Eta[P_Index[0]], P_LOOSE, PhoMode_);

                            ExcitedQuarkPlots[RunDependies[irun_]]->Fill(0.,weight*scaleFactor_);

                            for(int j=0;j<NJets;j++) 
                                ExcitedQuarkPlots[JetPtRunIdx[irun_]]->Fill(JetInfo.Pt[J_Index[j]],weight*scaleFactor_);
                        }

                        // control region
                        if(NPhotons==0&&NLeptons==2&&NLeptonsV==2&&NJets>=5){
                            float scaleFactor_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG))
                                scaleFactor_= 
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]],LepInfo.Eta[L_Index[0]], LepMode_)*
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[1]],LepInfo.Pt[L_Index[1]],LepInfo.Eta[L_Index[1]], LepMode_);
                            ExcitedQuarkPlots[_N2lep5jet+1+irun_]->Fill(0.,weight*scaleFactor_);
                        }
                        // signal region
                        if(NPhotons==2&&NLeptonsV==1){
                            float scaleFactor_ = 1.0;
                            float scaleFactorLep_ = 1.0;
                            float scaleFactorPho_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
                                scaleFactorLep_ = 
                                ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_);
                                scaleFactorPho_ = 
                                ScalePhotonID(PhotonInfo.Pt[P_Index[0]], PhotonInfo.Eta[P_Index[0]], P_LOOSE, PhoMode_)*
                                ScalePhotonID(PhotonInfo.Pt[P_Index[1]], PhotonInfo.Eta[P_Index[1]], P_LOOSE, PhoMode_);
                            }
                            scaleFactor_ = scaleFactorLep_ * scaleFactorPho_;

                            for(int j=0;j<NJets;j++) 
                                ExcitedQuarkPlots[JetPtRunIdx[irun_]+_JetPtRunAInSR-_JetPtRunA]->Fill(JetInfo.Pt[J_Index[j]],weight*scaleFactor_);
                        }
                    }
            }
        }// irun -- end --

        // For all runs
        for(int irun_=4;irun_<NRuns;irun_++){

            if(!strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
                if(!(EvtInfo.RunNo>=RunsRanges[irun_][0]&&EvtInfo.RunNo<=RunsRanges[irun_][1])) continue;
                if(!isGoodEvt(EvtInfo.RunNo,EvtInfo.LumiNo)) continue;

                // remove duplicate event
                map< pair<int, int> , int>::iterator evtitr;
                evtitr = evtlist[irun_].find( pair<int, int>(EvtInfo.RunNo, EvtInfo.EvtNo) );
                if( evtitr == evtlist[irun_].end() )
                    evtlist[irun_].insert( pair<pair<int, int>, int>(pair<int, int>(EvtInfo.RunNo, EvtInfo.EvtNo), 1));
                else
                    continue;
            }

            // Object selection -- start --
            int NMuons = 0;
            int NElectrons = 0;
            int NLeptons = 0;
            int NMuonsC = 0;
            int NElectronsC = 0;
            int NLeptonsC = 0;
            int NMuonsV = 0;
            int NElectronsV = 0;
            int NLeptonsV = 0;
            int NPhotons = 0;
            int NPhotonsNoISO = 0;
            int NJets = 0;
            int M_Index[MAX_LEPTONS/2];
            int E_Index[MAX_LEPTONS/2];
            int L_Index[MAX_LEPTONS/2*2];
            int MC_Index[MAX_LEPTONS/2];
            int EC_Index[MAX_LEPTONS/2];
            int LC_Index[MAX_LEPTONS/2*2];
            int MV_Index[MAX_LEPTONS/2];
            int EV_Index[MAX_LEPTONS/2];
            int LV_Index[MAX_LEPTONS/2*2];
            int P_Index[MAX_PHOTONS];
            int PnoISO_Index[MAX_PHOTONS];
            int J_Index[MAX_JETS];
            for(int m=0;m<MAX_LEPTONS/2;m++) M_Index[m] = 0;
            for(int e=0;e<MAX_LEPTONS/2;e++) E_Index[e] = 0;
            for(int l=0;l<MAX_LEPTONS/2*2;l++) L_Index[l] = 0;

            for(int m=0;m<MAX_LEPTONS/2;m++) MC_Index[m] = 0;
            for(int e=0;e<MAX_LEPTONS/2;e++) EC_Index[e] = 0;
            for(int l=0;l<MAX_LEPTONS/2*2;l++) LC_Index[l] = 0;
            for(int m=0;m<MAX_LEPTONS/2;m++) MV_Index[m] = 0;
            for(int e=0;e<MAX_LEPTONS/2;e++) EV_Index[e] = 0;
            for(int l=0;l<MAX_LEPTONS/2*2;l++) LV_Index[l] = 0;

            for(int p=0;p<MAX_PHOTONS;p++) P_Index[p] = 0;
            for(int p=0;p<MAX_PHOTONS;p++) PnoISO_Index[p] = 0;
            for(int j=0;j<MAX_JETS;j++) J_Index[j] = 0;
            RecoLooseLeptonSelection(EvtInfo,LepInfo,NMuonsV,MV_Index,NElectronsV,EV_Index,NLeptonsV,LV_Index);
            RecoLeptonSelection(EvtInfo,LepInfo,NMuons,M_Index,NElectrons,E_Index,NLeptons,L_Index);
            RecoLeptonForCleaningSelection(EvtInfo,LepInfo,NMuonsC,MC_Index,NElectronsC,EC_Index,NLeptonsC,LC_Index);

            UsingPhoISO = 10000000000.0;    // meaning no ISO reqirement
            RecoPhotonSelection(LepInfo,PhotonInfo,NMuonsC,MC_Index,NElectronsC,EC_Index,NPhotonsNoISO,PnoISO_Index,EvtInfo.RhoPU[0],P_LOOSE);
            UsingPhoISO = 1.0;
            RecoPhotonSelection(LepInfo,PhotonInfo,NMuonsC,MC_Index,NElectronsC,EC_Index,NPhotons,P_Index,EvtInfo.RhoPU[0],P_LOOSE);
            // there is a difference using standard lepton (muon) for photon cleaning
            //RecoPhotonSelection(LepInfo,PhotonInfo,NMuons,M_Index,NElectrons,E_Index,NPhotons,P_Index,EvtInfo.RhoPU[0],P_LOOSE);

            RecoJetSelection(LepInfo,JetInfo, PhotonInfo,NMuonsC,MC_Index,NElectronsC,EC_Index,NJets,J_Index,NPhotons,P_Index);

            // number of vertex -- start --
            int NOfVertex_=0;
            int V_Index[MAX_Vertices];
            for(int v=0;v<MAX_Vertices;v++) V_Index[v] = 0;
            for(int j=0;j<VertexInfo.Size;j++) {
                if(VertexInfo.isValid[j]!=1) continue; 
                if(VertexInfo.Type[j]!=0) continue; 
                if(VertexInfo.Ndof[j]<=4) continue; 
                if(fabs(VertexInfo.Rho[j])>=2) continue; 
                if(fabs(VertexInfo.z[j])>=24) continue; 
                V_Index[NOfVertex_] = j;
                NOfVertex_++;
            }
            if(NOfVertex_<1) continue;

            // 5) MM
            if((NLeptons>=1&&NJets>=4)) {
                if(JetInfo.Pt[J_Index[0]]>JetPtOrderCuts[irun_][0]) 
                if(JetInfo.Pt[J_Index[1]]>JetPtOrderCuts[irun_][1]) 
                    if(JetInfo.Pt[J_Index[2]]>JetPtOrderCuts[irun_][2]){ 
                        bool IsMuon =0;
                        if (LepInfo.LeptonType[L_Index[0]]==13) IsMuon =1;
                        if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG)) {
                            ExcitedQuarkPlots[_NOfVertexb4PU1l4j]->Fill(NOfVertex_,weight_*LUMINOSITYForRuns[irun_]/LUMINOSITY);

                            float WeightCorrectionOnTstar = 1.0;
                            if(sample_idx>=Tgamma500GeV&&sample_idx<=TGluonTgamma1500GeV){
                                int MCTruthChannel_ = MCTruthChannel(GenInfo);
                                WeightCorrectionOnTstar = CorrectionOnTstar(sample_idx,MCTruthChannel_);
                            }

                            if(NPhotons>=1&&NLeptonsV==1){
                                ExcitedQuarkPlots[_NOfVertexb4PU1p1l4j]->Fill(NOfVertex_,
                                        weight_*LUMINOSITYForRuns[irun_]/LUMINOSITY);
                                ExcitedQuarkPlots[_NOfVertexb4PU1p1l4jWithTstarCorrection]->Fill(NOfVertex_,
                                        WeightCorrectionOnTstar*weight_*LUMINOSITYForRuns[irun_]/LUMINOSITY);
                            }

                            if((EvtInfo.TrueIT[0])>=60){
                                std::cout<<"[PU WARNING] N(vertex = "<<EvtInfo.TrueIT[0]<<") >= 60"<<std::endl;
                            }else{
                                weight = WeightCorrectionOnTstar*weight_*PileUpReweight[irun_]->GetBinContent(PileUpReweight[irun_]->GetXaxis()->FindBin(EvtInfo.TrueIT[0]))*LUMINOSITYForRuns[irun_]/LUMINOSITY;
                            }


                            ExcitedQuarkPlots[_NOfVertexafPU1l4j]->Fill(NOfVertex_,weight);
                            if(NPhotons>=1&&NLeptonsV==1)
                                ExcitedQuarkPlots[_NOfVertexafPU1p1l4j]->Fill(NOfVertex_,weight);
                            // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiTopRefEventSel#Trigger_Efficiencies
                            std::vector<double> WeightTopTriggerEfficiency = weight_provider[irun_]->get_weight(
                                    LepInfo.Pt[L_Index[0]],LepInfo.Eta[L_Index[0]],
                                    JetInfo.Pt[J_Index[3]],JetInfo.Eta[J_Index[3]],
                                    NOfVertex_,
                                    NJets,
                                    IsMuon,
                                    TopTriggerEfficiencyProvider::NOMINAL);
                            double uncOnTrig = 0.;
                            if(RunStatus_ == UncTrigPlus) {
                                uncOnTrig = 1.0;
                                //uncOnTrig = 0.01;
                            }else if (RunStatus_ == UncTrigMinus){
                                uncOnTrig = -1.0;
                                //uncOnTrig = -0.01;
                            }
                            weight *= (WeightTopTriggerEfficiency[0] + uncOnTrig*WeightTopTriggerEfficiency[1]);
                            //weight *= (WeightTopTriggerEfficiency[0] + uncOnTrig);
                            ExcitedQuarkPlots[_NOfVertexafPUTrgeff1l4j]->Fill(NOfVertex_,weight);
                            if(NLeptonsV==1){
                                float scaleFactor_ =  
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_);

                                float TopPtReweighting = 1.0;
                                if (sample_idx == 
                                TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1){
                                    int idxes[2] = {-1,-1};
                                    for(int iTop=0;iTop<GenInfo.Size;iTop++){
                                        if(GenInfo.PdgID[iTop] == 6) idxes[0] = iTop; 
                                        if(GenInfo.PdgID[iTop] == -6) idxes[1] = iTop; 
                                    }   
                                    if(idxes[0]!=-1&&idxes[1]!=-1){
                                        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
                                        TopPtReweighting = 
                                        sqrt(exp(0.159-0.00141*GenInfo.Pt[idxes[0]])*exp(0.159-0.00141*GenInfo.Pt[idxes[1]]));
                                        TopPtReweighting = pow(TopPtReweighting,TopPtMode_);
                                    }else{
                                        std::cout<<"[WARNING] no Gen top"<<std::endl;
                                    }

                                }

                                if (sample_idx >= 
                                        TTJets_FullLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v2 && 
                                        (sample_idx<QCD_Pt_15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2||
                                         sample_idx>QCD_Pt_800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1
                                        )
                                        ){
                                    double weightTemp = weight*scaleFactor_*TopPtReweighting ;

                                    // Denominator
                                    for(int p=0;p<NPhotonsNoISO;p++){
                                        int objIdx = PnoISO_Index[p];
                                        int mctruth_ = ObjectMCTruthForRatio(objIdx, GenInfo, PhotonInfo);
                                        int cat_tag = -1;
                                        if(mctruth_==0){
                                            cat_tag = 0;    // prompt photon
                                        }else if(mctruth_==6||mctruth_==12||mctruth_==18){
                                            cat_tag = 1;    // decay in flight/ISR/FSR
                                        }else if(mctruth_==15){
                                            cat_tag = 2;    // gluon
                                        }else if(mctruth_==5){
                                            cat_tag = 3;    // light quark
                                        }else if(mctruth_==3 || mctruth_==4){
                                            cat_tag = 4;    // heavy quark
                                        }else if(mctruth_==9){
                                            cat_tag = 5;    // electron
                                        }else if(mctruth_==8){
                                            cat_tag = 6;    // muon
                                        }else if(mctruth_==7){
                                            cat_tag = 7;    // tau
                                        }
                                        if(cat_tag==-1) continue;
                                        hMCBGpTandEtaDenominator[cat_tag]->Fill(PhotonInfo.Pt[objIdx], 
                                                PhotonInfo.Eta[objIdx], weightTemp);
                                        hMCBGpTDenominator[cat_tag]->Fill(PhotonInfo.Pt[objIdx], weightTemp);
                                        hMCBGEtaDenominator[cat_tag]->Fill(PhotonInfo.Eta[objIdx], weightTemp);
                                    }

                                    // Numerator
                                    for(int p=0;p<NPhotons;p++){
                                        int objIdx = P_Index[p];
                                        int mctruth_ = ObjectMCTruthForRatio(objIdx, GenInfo, PhotonInfo);
                                        int cat_tag = -1;
                                        if(mctruth_==0){
                                            cat_tag = 0;    // photon
                                        }else if(mctruth_==6||mctruth_==12||mctruth_==18){
                                            cat_tag = 1;    // decay in flight/ISR/FSR
                                        }else if(mctruth_==15){
                                            cat_tag = 2;    // gluon
                                        }else if(mctruth_==5){
                                            cat_tag = 3;    // light quark
                                        }else if(mctruth_==3 || mctruth_==4){
                                            cat_tag = 4;    // heavy quark
                                        }else if(mctruth_==9){
                                            cat_tag = 5;    // electron
                                        }else if(mctruth_==8){
                                            cat_tag = 6;    // muon
                                        }else if(mctruth_==7){
                                            cat_tag = 7;    // tau
                                        }
                                        if(cat_tag==-1) continue;
                                        hMCBGpTandEtaNumerator[cat_tag]->Fill(PhotonInfo.Pt[objIdx], 
                                                PhotonInfo.Eta[objIdx], weightTemp);
                                        hMCBGpTNumerator[cat_tag]->Fill(PhotonInfo.Pt[objIdx], weightTemp);
                                        hMCBGEtaNumerator[cat_tag]->Fill(PhotonInfo.Eta[objIdx], weightTemp);
                                    }
                                }
                            }
                            if(NPhotons>=1&&NLeptonsV==1){
                                ExcitedQuarkPlots[_NOfVertexafPUTrgeff1p1l4j]->Fill(NOfVertex_,weight);
                                float scaleFactor_ =  
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_);
                                ExcitedQuarkPlots[_NOfVertexafPUTrgeff1p1l4jafLS]->Fill(NOfVertex_,
                                        weight*scaleFactor_ );
                                scaleFactor_= 
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_)*
                                    ScalePhotonID(PhotonInfo.Pt[P_Index[0]], PhotonInfo.Eta[P_Index[0]], P_LOOSE, PhoMode_);
                                ExcitedQuarkPlots[_NOfVertexafPUTrgeff1p1l4jafLSPS]->Fill(NOfVertex_,
                                        weight*scaleFactor_ );

                                float TopPtReweighting = 1.0;
                                if (sample_idx == 
                                TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1){
                                    int idxes[2] = {-1,-1};
                                    for(int iTop=0;iTop<GenInfo.Size;iTop++){
                                        if(GenInfo.PdgID[iTop] == 6) idxes[0] = iTop; 
                                        if(GenInfo.PdgID[iTop] == -6) idxes[1] = iTop; 
                                    }   
                                    if(idxes[0]!=-1&&idxes[1]!=-1){
                                        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
                                        TopPtReweighting = 
                                        sqrt(exp(0.159-0.00141*GenInfo.Pt[idxes[0]])*exp(0.159-0.00141*GenInfo.Pt[idxes[1]]));
                                        TopPtReweighting = pow(TopPtReweighting,TopPtMode_);
                                    }else{
                                        std::cout<<"[WARNING] no Gen top"<<std::endl;
                                    }

                                }
                                ExcitedQuarkPlots[_NOfVertexafPUTrgeff1p1l4jafLSPSTop]->Fill(NOfVertex_,
                                        weight*scaleFactor_*TopPtReweighting );
                                weight = weight*TopPtReweighting;
                            }
                        }else{
                            ExcitedQuarkPlots[_NOfVertexb4PU1l4j]->Fill(NOfVertex_);
                            ExcitedQuarkPlots[_NOfVertexafPU1l4j]->Fill(NOfVertex_);
                            ExcitedQuarkPlots[_NOfVertexafPUTrgeff1l4j]->Fill(NOfVertex_);
                            if(NPhotons>=1&&NLeptonsV==1){
                                ExcitedQuarkPlots[_NOfVertexb4PU1p1l4j]->Fill(NOfVertex_);
                                ExcitedQuarkPlots[_NOfVertexafPU1p1l4j]->Fill(NOfVertex_);
                                ExcitedQuarkPlots[_NOfVertexafPUTrgeff1p1l4j]->Fill(NOfVertex_);
                                ExcitedQuarkPlots[_NOfVertexafPUTrgeff1p1l4jafLS]->Fill(NOfVertex_);
                                ExcitedQuarkPlots[_NOfVertexafPUTrgeff1p1l4jafLSPS]->Fill(NOfVertex_);
                                ExcitedQuarkPlots[_NOfVertexafPUTrgeff1p1l4jafLSPSTop]->Fill(NOfVertex_);
                            }
                        }

                        if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
                            ExcitedQuark2DPlots[_NjetNlepMC]->Fill(NJets,NLeptons,weight);
                        }else{
                            ExcitedQuark2DPlots[_NjetNlepData]->Fill(NJets,NLeptons,weight);
                        }

                        // loose regions
                        ExcitedQuarkPlots[_N1lep4jet]->Fill(0.,weight);
                        if(NPhotons>=1&&NLeptonsV==1){
                            float scaleFactor_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG))
                                scaleFactor_= 
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_)*
                                    ScalePhotonID(PhotonInfo.Pt[P_Index[0]], PhotonInfo.Eta[P_Index[0]], P_LOOSE, PhoMode_);

                            ExcitedQuarkPlots[_N1p1lep4jet]->Fill(0.,weight*scaleFactor_);

                            for(int j=0;j<NJets;j++) ExcitedQuarkPlots[_JetPt]->Fill(JetInfo.Pt[J_Index[j]],weight*scaleFactor_);
                            for(int j=0;j<4;j++) ExcitedQuarkPlots[_Jet1stPt+j]->Fill(JetInfo.Pt[J_Index[j]],weight*scaleFactor_);
                            for(int p=0;p<NPhotons;p++) 
                                ExcitedQuarkPlots[_PhotonPt]->Fill(PhotonInfo.Pt[P_Index[p]],weight*scaleFactor_);
                            for(int e=0;e<NElectrons;e++) 
                                ExcitedQuarkPlots[_ElectronPt]->Fill(LepInfo.Pt[E_Index[e]],weight*scaleFactor_);
                            for(int m=0;m<NMuons;m++) 
                                ExcitedQuarkPlots[_MuonPt]->Fill(LepInfo.Pt[M_Index[m]],weight*scaleFactor_);

                            double _ST_ = 0.;
                            for(int j=0;j<NJets;j++) _ST_ += JetInfo.Pt[J_Index[j]];
                            for(int m=0;m<NMuons;m++) _ST_ += LepInfo.Pt[M_Index[m]];
                            for(int e=0;e<NElectrons;e++) _ST_ += LepInfo.Pt[E_Index[e]];
                            for(int p=0;p<NPhotons;p++) _ST_ += PhotonInfo.Pt[P_Index[p]];
                            _ST_+=EvtInfo.PFMET;
                            ExcitedQuarkPlots[_ST]->Fill(_ST_,weight*scaleFactor_);
                            ExcitedQuarkPlots[_Nj]->Fill(NJets,weight*scaleFactor_);
                            ExcitedQuarkPlots[_Np]->Fill(NPhotons,weight*scaleFactor_);
                            ExcitedQuarkPlots[_Ne]->Fill(NElectrons,weight*scaleFactor_);
                            ExcitedQuarkPlots[_Nu]->Fill(NMuons,weight*scaleFactor_);
                            ExcitedQuarkPlots[_MET]->Fill(EvtInfo.PFMET,weight*scaleFactor_);

                            fprintf(pdffileCR,"%i %i %f %f %f %f %f %f\n",
                                    EvtInfo.PDFid1,
                                    EvtInfo.PDFid2,
                                    EvtInfo.PDFx1,
                                    EvtInfo.PDFx2,
                                    EvtInfo.PDFscale,
                                    EvtInfo.PDFv1,
                                    EvtInfo.PDFv2,
                                    weight*scaleFactor_
                                    );

                            int MCTruthChannel_ = MCTruthChannel(GenInfo);
                            if(MCTruthChannel_!=-1)
                                MCTruthChannels[MCTruthChannel_]++;
                            MCTruthChannels[3]++;
                        }

                        Nr_template=-1;        
                        Nj_template=-1;
                        Nl_template=-1;
                        Nlloose_template=-1;
                        Nrloose_template=-1;
                        weight_template = -1;
                        Mass_template=-1;
                        category_template=-1;

                        Nr_template = NPhotons;
                        Nj_template = NJets;        
                        Nl_template = NLeptons;
                        Nlloose_template = NLeptonsV;
                        Nrloose_template = -1;  
                        MCTruth_tmp[0] = -1;
                        MCTruth_tmp[1] = -1;
                        MCTruth_tmp[2] = -1;
                        MCTruth1 = -1;
                        MCTruth2 = -1;
                        MCTruth3 = -1;

                        // control region
                        if(NPhotons==1&&NLeptons==1&&NLeptonsV==1&&NJets>=5){
                            float scaleFactor_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG))
                                scaleFactor_= 
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_)*
                                    ScalePhotonID(PhotonInfo.Pt[P_Index[0]], PhotonInfo.Eta[P_Index[0]], P_LOOSE, PhoMode_);
                            ExcitedQuarkPlots[_N1p1lep5jet]->Fill(0.,weight*scaleFactor_);
                            if(RunForCR_Shape){
                                if(sample_idx<TTJets_FullLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v2||
                                        sample_idx>ZH_HToZZTo4L_M_126_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7C_v1){
                                    category_template = 0; 
                                    Mass_template =                 
                                        RecoFakeTstar(EvtInfo,GenInfo,LepInfo,JetInfo,PhotonInfo,NLeptons,L_Index,NJets,  
                                                J_Index,NPhotons,  P_Index,MCTruth_tmp,category_template);
                                    ExcitedQuarkPlots[_MtstarChi2_1r1l5j]->Fill( Mass_template ,weight*scaleFactor_);

                                    weight_template = weight*scaleFactor_;
                                    MCTruth1 = MCTruth_tmp[0];
                                    MCTruth2 = MCTruth_tmp[1];
                                    MCTruth3 = MCTruth_tmp[2];
                                }
                            }
                            Nr_template = NPhotons;
                            Nj_template = NJets;        
                            Nl_template = NLeptons;
                            MCTemplatesTree->Fill();
                        }
                        if(NPhotons==0&&NLeptons==1&&NLeptonsV==1&&NJets>=6){
                            float scaleFactor_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG))
                                scaleFactor_= 
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_);
                            ExcitedQuarkPlots[_N1lep6jet]->Fill(0.,weight*scaleFactor_);
                            if(RunForCR_Shape){
                                if(sample_idx<TTJets_FullLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v2||
                                        sample_idx>ZH_HToZZTo4L_M_126_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7C_v1){
                                    category_template = 1; 
                                    Mass_template =
                                        RecoFakeTstar(EvtInfo,GenInfo,LepInfo,JetInfo,PhotonInfo,NLeptons,L_Index,NJets,  
                                                J_Index,NPhotons,  P_Index,MCTruth_tmp,category_template);
                                    ExcitedQuarkPlots[_MtstarChi2_0r1l6j]->Fill(Mass_template,weight*scaleFactor_);
                                    weight_template = weight*scaleFactor_;
                                    MCTruth1 = MCTruth_tmp[0];
                                    MCTruth2 = MCTruth_tmp[1];
                                    MCTruth3 = MCTruth_tmp[2];
                                }
                            }
                            Nr_template = NPhotons;
                            Nj_template = NJets;        
                            Nl_template = NLeptons;
                            MCTemplatesTree->Fill();
                        }
                        if(NPhotons==1&&NLeptons==2&&NLeptonsV==2&&NJets>=4){
                            float scaleFactor_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG))
                                scaleFactor_= 
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_)*
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[1]],LepInfo.Pt[L_Index[1]], LepInfo.Eta[L_Index[1]], LepMode_)*
                                    ScalePhotonID(PhotonInfo.Pt[P_Index[0]], PhotonInfo.Eta[P_Index[0]], P_LOOSE, PhoMode_);
                            ExcitedQuarkPlots[_N1p2lep4jet]->Fill(0.,weight*scaleFactor_);
                            if(RunForCR_Shape){
                                if(sample_idx<TTJets_FullLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v2||
                                        sample_idx>ZH_HToZZTo4L_M_126_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7C_v1){
                                    category_template = 2;
                                    Mass_template =
                                        RecoFakeTstar(EvtInfo,GenInfo,LepInfo,JetInfo,PhotonInfo,NLeptons,L_Index,NJets,  
                                                J_Index,NPhotons,  P_Index,MCTruth_tmp,category_template);
                                    ExcitedQuarkPlots[_MtstarChi2_1r2l4j]->Fill(Mass_template,weight*scaleFactor_);
                                    weight_template = weight*scaleFactor_;
                                    MCTruth1 = MCTruth_tmp[0];
                                    MCTruth2 = MCTruth_tmp[1];
                                    MCTruth3 = MCTruth_tmp[2];
                                }
                            }
                            Nr_template = NPhotons;
                            Nj_template = NJets;        
                            Nl_template = NLeptons;
                            MCTemplatesTree->Fill();
                        }
                        if(NPhotons==0&&NLeptons>=3&&NLeptonsV==3&&NJets>=4){
                            float scaleFactor_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG))
                                scaleFactor_= 
                                ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_)*
                                ScaleLeptonID(LepInfo.LeptonType[L_Index[1]],LepInfo.Pt[L_Index[1]], LepInfo.Eta[L_Index[1]], LepMode_)*
                                ScaleLeptonID(LepInfo.LeptonType[L_Index[2]],LepInfo.Pt[L_Index[2]], LepInfo.Eta[L_Index[2]], LepMode_);
                            ExcitedQuarkPlots[_N3lep4jet]->Fill(0.,weight*scaleFactor_);
                            if(RunForCR_Shape){
                                if(sample_idx<TTJets_FullLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v2||
                                        sample_idx>ZH_HToZZTo4L_M_126_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7C_v1){
                                    category_template = 3;
                                    Mass_template =
                                        RecoFakeTstar(EvtInfo,GenInfo,LepInfo,JetInfo,PhotonInfo,NLeptons,L_Index,NJets,  
                                                J_Index,NPhotons,  P_Index,MCTruth_tmp,category_template);
                                    ExcitedQuarkPlots[_MtstarChi2_0r3l4j]->Fill(Mass_template,weight*scaleFactor_);
                                    weight_template = weight*scaleFactor_;
                                    MCTruth1 = MCTruth_tmp[0];
                                    MCTruth2 = MCTruth_tmp[1];
                                    MCTruth3 = MCTruth_tmp[2];
                                }
                            }
                            Nr_template = NPhotons;
                            Nj_template = NJets;        
                            Nl_template = NLeptons;
                            MCTemplatesTree->Fill();
                        }
                        if(NPhotons==0&&NLeptons==2&&NLeptonsV==2&&NJets>=5){
                            float scaleFactor_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG))
                                scaleFactor_= 
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]],LepInfo.Eta[L_Index[0]], LepMode_)*
                                    ScaleLeptonID(LepInfo.LeptonType[L_Index[1]],LepInfo.Pt[L_Index[1]],LepInfo.Eta[L_Index[1]], LepMode_);
                            ExcitedQuarkPlots[_N2lep5jet]->Fill(0.,weight*scaleFactor_);
                            if(RunForCR_Shape){
                                if(sample_idx<TTJets_FullLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v2||
                                        sample_idx>ZH_HToZZTo4L_M_126_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7C_v1){
                                    category_template = 4;
                                    Mass_template =
                                        RecoFakeTstar(EvtInfo,GenInfo,LepInfo,JetInfo,PhotonInfo,NLeptons,L_Index,NJets,  
                                                J_Index,NPhotons,  P_Index,MCTruth_tmp,category_template);
                                    ExcitedQuarkPlots[_MtstarChi2_0r2l5j]->Fill(Mass_template,weight*scaleFactor_);
                                    weight_template = weight*scaleFactor_;
                                    MCTruth1 = MCTruth_tmp[0];
                                    MCTruth2 = MCTruth_tmp[1];
                                    MCTruth3 = MCTruth_tmp[2];
                                }
                            }

                            Nr_template = NPhotons;
                            Nj_template = NJets;        
                            Nl_template = NLeptons;
                            MCTemplatesTree->Fill();
                        }
                        // signal region
                        //if(NPhotons>=2&&NLeptonsV==1){
                        if(NPhotons==2&&NLeptonsV==1){
                            float scaleFactor_ = 1.0;
                            float scaleFactorLep_ = 1.0;
                            float scaleFactorPho_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
                                scaleFactorLep_ = 
                                ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_);
                                scaleFactorPho_ = 
                                ScalePhotonID(PhotonInfo.Pt[P_Index[0]], PhotonInfo.Eta[P_Index[0]], P_LOOSE, PhoMode_)*
                                ScalePhotonID(PhotonInfo.Pt[P_Index[1]], PhotonInfo.Eta[P_Index[1]], P_LOOSE, PhoMode_);
                            }
                            scaleFactor_ = scaleFactorLep_ * scaleFactorPho_;
                            ExcitedQuarkPlots[_N2p1lep4jet]->Fill(0.,weight*scaleFactor_);
                            double *Mtstar_ =  RecoTstar(
                                    EvtInfo,
                                    GenInfo,
                                    LepInfo,
                                    JetInfo,
                                    PhotonInfo,
                                    NLeptons,L_Index,
                                    NJets,  J_Index,
                                    NPhotons,  P_Index,MCTruth_tmp, 0);
                            Mass_template = Mtstar_[1];
                            ExcitedQuarkPlots[_MtstarChi2_2r1l4jnoLSPS]->Fill(Mass_template,weight);
                            ExcitedQuarkPlots[_MtstarChi2_2r1l4jnoPS]->Fill(Mass_template,weight*scaleFactorLep_);
                            ExcitedQuarkPlots[_MtstarChi2_2r1l4j]->Fill(Mass_template,weight*scaleFactor_);
                            if(!strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
                                ExcitedQuarkPlots[_MtstarChi2_2r1l4jReweight1]->Fill(Mass_template,weight*scaleFactor_);
                                ExcitedQuarkPlots[_MtstarChi2_2r1l4jReweight2]->Fill(Mass_template,weight*scaleFactor_);
                            }
                            category_template = 10;

                            weight_template = weight*scaleFactor_;
                            MCTruth1 = MCTruth_tmp[0];
                            MCTruth2 = MCTruth_tmp[1];
                            MCTruth3 = MCTruth_tmp[2];

                            for(int j=0;j<NJets;j++) 
                                ExcitedQuarkPlots[_JetPtInSR]->Fill(JetInfo.Pt[J_Index[j]],weight*scaleFactor_);
                            for(int j=0;j<4;j++) 
                                ExcitedQuarkPlots[_Jet1stPtInSR+j]->Fill(JetInfo.Pt[J_Index[j]],weight*scaleFactor_);
                            for(int p=0;p<NPhotons;p++) 
                                ExcitedQuarkPlots[_PhotonPtInSR]->Fill(PhotonInfo.Pt[P_Index[p]],weight*scaleFactor_);
                            for(int e=0;e<NElectrons;e++) 
                                ExcitedQuarkPlots[_ElectronPtInSR]->Fill(LepInfo.Pt[E_Index[e]],weight*scaleFactor_);
                            for(int m=0;m<NMuons;m++) 
                                ExcitedQuarkPlots[_MuonPtInSR]->Fill(LepInfo.Pt[M_Index[m]],weight*scaleFactor_);

                            double _ST_ = 0.;
                            for(int j=0;j<NJets;j++) _ST_ += JetInfo.Pt[J_Index[j]];
                            for(int m=0;m<NMuons;m++) _ST_ += LepInfo.Pt[M_Index[m]];
                            for(int e=0;e<NElectrons;e++) _ST_ += LepInfo.Pt[E_Index[e]];
                            for(int p=0;p<NPhotons;p++) _ST_ += PhotonInfo.Pt[P_Index[p]];
                            _ST_+=EvtInfo.PFMET;
                            ExcitedQuarkPlots[_STInSR]->Fill(_ST_,weight*scaleFactor_);
                            ExcitedQuarkPlots[_NjInSR]->Fill(NJets,weight*scaleFactor_);
                            ExcitedQuarkPlots[_NpInSR]->Fill(NPhotons,weight*scaleFactor_);
                            ExcitedQuarkPlots[_NeInSR]->Fill(NElectrons,weight*scaleFactor_);
                            ExcitedQuarkPlots[_NuInSR]->Fill(NMuons,weight*scaleFactor_);
                            ExcitedQuarkPlots[_METInSR]->Fill(EvtInfo.PFMET,weight*scaleFactor_);

                            fprintf(pdffileSR,"%i %i %f %f %f %f %f %f\n",
                                    EvtInfo.PDFid1,
                                    EvtInfo.PDFid2,
                                    EvtInfo.PDFx1,
                                    EvtInfo.PDFx2,
                                    EvtInfo.PDFscale,
                                    EvtInfo.PDFv1,
                                    EvtInfo.PDFv2,
                                    weight*scaleFactor_
                                    );

                            if(!strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
                                std::cout<<"[Selected "<<EvtInfo.RunNo <<", "<<EvtInfo.LumiNo<<", "
                                    <<EvtInfo.EvtNo <<"] w/ (Lep, nj, hmass, lmass, ST)=("<<
                                    LepInfo.LeptonType[L_Index[0]] <<", "<<NJets<<", "<<Mtstar_[0]<<", "<<Mtstar_[1]<<", "<<
                                    _ST_<<")"<<std::endl;
                                for(int j=0;j<NJets;j++) 
                                    std::cout<<"        [Selected jet#"<<j<<"] w/ (pT, eta, phi)=("<<
                                        JetInfo.Pt[J_Index[j]]<<", "<<JetInfo.Eta[J_Index[j]]<<", "<<
                                        JetInfo.Phi[J_Index[j]]<<")"<<std::endl;
                                for(int l=0;l<NLeptons;l++)
                                    std::cout<<"        [Selected lep#"<<l<<"] w/ (pT, eta, phi)=("<<
                                        LepInfo.Pt[L_Index[l]]<<", "<<LepInfo.Eta[L_Index[l]]<<", "<<
                                        LepInfo.Phi[L_Index[l]]<<")"<<std::endl;
                                for(int p=0;p<NPhotons;p++) 
                                    std::cout<<"        [Selected pho#"<<p<<"] w/ (pT, eta, phi)=("<<
                                        PhotonInfo.Pt[P_Index[p]]<<", "<<PhotonInfo.Eta[P_Index[p]]<<", "<<
                                        PhotonInfo.Phi[P_Index[p]]<<")"<<std::endl;
                                std::cout<<"        [Selected MET] w/ (pT, phi)=("<<
                                    EvtInfo.PFMET<<", "<<EvtInfo.PFMETPhi<<")"<<std::endl;                      
                            }
                            Nr_template = NPhotons;
                            Nj_template = NJets;        
                            Nl_template = NLeptons;
                            MCTemplatesTree->Fill();
                        }
                        // Loose region x ISO ratio --> mimic signal region
                        // [be careful] This region will interfere with other region, and then reset is needed.
                        Nr_template=-1;        
                        Nj_template=-1;
                        Nl_template=-1;
                        Nlloose_template=-1;
                        Nrloose_template=-1;
                        weight_template = -1;
                        Mass_template=-1;
                        category_template=-1;

                        Nr_template = NPhotons;
                        Nj_template = NJets;        
                        Nl_template = NLeptons;
                        Nlloose_template = NLeptonsV;
                        Nrloose_template = -1;  
                        MCTruth_tmp[0] = -1;
                        MCTruth_tmp[1] = -1;
                        MCTruth_tmp[2] = -1;
                        MCTruth1 = -1;
                        MCTruth2 = -1;
                        MCTruth3 = -1;
                        if(NPhotonsNoISO==2&&NLeptonsV==1&&strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
                            float scaleFactor_ = 1.0;
                            float scaleFactorLep_ = 1.0;
                            float scaleFactorPho_ = 1.0;
                            if(strcmp(SAMPLE[sample_idx].tag,DATA_TAG)){
                                scaleFactorLep_ = 
                                ScaleLeptonID(LepInfo.LeptonType[L_Index[0]],LepInfo.Pt[L_Index[0]], LepInfo.Eta[L_Index[0]], LepMode_);
                                scaleFactorPho_ = 
                                    ScalePhotonID(PhotonInfo.Pt[PnoISO_Index[0]], PhotonInfo.Eta[PnoISO_Index[0]], P_LOOSE, PhoMode_)*
                                    ScalePhotonID(PhotonInfo.Pt[PnoISO_Index[1]], PhotonInfo.Eta[PnoISO_Index[1]], P_LOOSE, PhoMode_);
                            }
                            scaleFactor_ = scaleFactorLep_ * scaleFactorPho_;
                            double *Mtstar_ =  RecoTstar(
                                    EvtInfo,
                                    GenInfo,
                                    LepInfo,
                                    JetInfo,
                                    PhotonInfo,
                                    NLeptons,L_Index,
                                    NJets,  J_Index,
                                    NPhotonsNoISO,  PnoISO_Index,MCTruth_tmp, 0);
                            Mass_template = Mtstar_[1];
                            category_template = 5;

                            MCTruth1 = MCTruth_tmp[0];
                            MCTruth2 = MCTruth_tmp[1];
                            MCTruth3 = MCTruth_tmp[2];

                            // Find MC truth for 2 photons
                            int CAT_TAG[NPhotonsNoISO];
                            int BinsFor2Ps[NPhotonsNoISO];
                            for(int p=0;p<NPhotonsNoISO;p++){
                                CAT_TAG[p] = -1;
                                BinsFor2Ps[p] = -1;
                                int objIdx = PnoISO_Index[p];
                                int mctruth_ = ObjectMCTruthForRatio(objIdx, GenInfo, PhotonInfo);
                                int cat_tag = -1;
                                if(mctruth_==0){
                                    cat_tag = 0;    // prompt photon
                                }else if(mctruth_==6||mctruth_==12||mctruth_==18){
                                    cat_tag = 1;    // decay in flight/ISR/FSR
                                }else if(mctruth_==15){
                                    cat_tag = 2;    // gluon
                                }else if(mctruth_==5){
                                    cat_tag = 3;    // light quark
                                }else if(mctruth_==3 || mctruth_==4){
                                    cat_tag = 4;    // heavy quark
                                }else if(mctruth_==9||mctruth_==8){
                                    cat_tag = 5;    // electron or muon
                                //}else if(mctruth_==8){
                                //    cat_tag = 6;    // muon
                                }else if(mctruth_==7){
                                    cat_tag = 7;    // tau
                                }
                                CAT_TAG[p] = cat_tag;
                                BinsFor2Ps[p] = 
                                    hRatioMCBGpTandEta[0]->FindBin(PhotonInfo.Pt[objIdx],PhotonInfo.Eta[objIdx]);
                            }
                            if(CAT_TAG[0]!=-1&&CAT_TAG[1]!=-1){
                                weight_template = weight*scaleFactor_*
                                    hRatioMCBGpTandEta[CAT_TAG[0]]->GetBinContent(BinsFor2Ps[0])*
                                    hRatioMCBGpTandEta[CAT_TAG[1]]->GetBinContent(BinsFor2Ps[1]);
                            }else{
                                weight_template = -1;
                            }
                            if(weight_template != -1)
                                ExcitedQuarkPlots[_MtstarChi2_2r1l4jReweight1]->Fill(Mass_template,weight_template);

                            if(weight_template == -1)
                                ExcitedQuarkPlots[_MtstarChi2_2r1l4jReweight2]->Fill(Mass_template,weight*scaleFactor_);
                            else
                                ExcitedQuarkPlots[_MtstarChi2_2r1l4jReweight2]->Fill(Mass_template,weight_template);

                            Nr_template = NPhotonsNoISO;
                            Nj_template = NJets;        
                            Nl_template = NLeptons;
                            MCTemplatesTree->Fill();
                        }
                    }
            }
        }// irun for All Runs -- end --
    }
    MCTemplatesTree->Write();
    root->Delete();
    fclose(pdffileCR);
    fclose(pdffileSR);
    //delete weight_provider;
}

void integralWithCut(TH1F *h1,double cut){
	double leftNumber = 0;
	double rightNumber = 0;

	for(int i=0;i<h1->GetXaxis()->GetNbins();i++){
		if(h1->GetXaxis()->GetBinCenter(i)<cut){
			leftNumber += h1->GetBinContent(i);
		}else{
			rightNumber += h1->GetBinContent(i);
		}
	}
	std::cout<<"cut ( "<<cut<<" ) , rightNumber : "<<rightNumber<<" , leftNumber : "<<leftNumber<<std::endl;
}

void proj_anyregions_shape_loose_norunCR(){

    //if(RunStatus_==Normal){
    if(true){   // also consider uncertainty on the control samples
        RunForCR_Shape = true;
        RunForMoreStat = true;
    }

	double time1,time2,time3; time1=time(NULL);
	char buffer_1[128];

	// Get ratio for background estimation
	TFile *RatioFile = new TFile("interface/Ratio.root");
	hRatioData = (TH2D*) RatioFile->Get("hData_1tp1lep4jet");
	hRatioMCBG = (TH2D*) RatioFile->Get("hMCBG_1tp1lep4jet");


    // Get Efficiency(ISO) for weighting M(t+pho) without ISO
    for(int i=0;i<Nratio;i++){// [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
        sprintf(buffer_1,"hMCBGpTandEtaNumerator_%s",RatioNames[i].c_str());
        hMCBGpTandEtaNumerator[i] = new TH2D(buffer_1,"",400,0,400,100,-2.5,2.5);    
        hMCBGpTandEtaNumerator[i]->Sumw2();

        sprintf(buffer_1,"hMCBGpTandEtaDenominator_%s",RatioNames[i].c_str());
        hMCBGpTandEtaDenominator[i] = new TH2D(buffer_1,"",400,0,400,100,-2.5,2.5);  
        hMCBGpTandEtaDenominator[i]->Sumw2();

        sprintf(buffer_1,"hMCBGpTNumerator_%s",RatioNames[i].c_str());
        hMCBGpTNumerator[i] = new TH1D(buffer_1,"",400,0,400); 
        hMCBGpTNumerator[i]->Sumw2();

        sprintf(buffer_1,"hMCBGpTDenominator_%s",RatioNames[i].c_str());
        hMCBGpTDenominator[i] = new TH1D(buffer_1,"",400,0,400);    
        hMCBGpTDenominator[i]->Sumw2();
        
        sprintf(buffer_1,"hMCBGEtaNumerator_%s",RatioNames[i].c_str());
        hMCBGEtaNumerator[i] = new TH1D(buffer_1,"",100,-2.5,2.5);    
        hMCBGEtaNumerator[i]->Sumw2();

        sprintf(buffer_1,"hMCBGEtaDenominator_%s",RatioNames[i].c_str());
        hMCBGEtaDenominator[i] = new TH1D(buffer_1,"",100,-2.5,2.5);    
        hMCBGEtaDenominator[i]->Sumw2();
    }
    // Get ISO ratio
	TFile *ISORatioFile = new TFile("interface/ISO_ratio.root");
    for(int i=0;i<Nratio;i++){// [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
        sprintf(buffer_1,"hRatioMCBGpTandEta_%s",RatioNames[i].c_str());
        hRatioMCBGpTandEta[i] = (TH2D*) ISORatioFile->Get(buffer_1);
    }

    // Trigger efficiency
    for(int i=0;i<NRuns;i++){
        weight_provider[i] = new TopTriggerEfficiencyProvider();
        if(i==0) weight_provider[i]->setLumi(TopTriggerEfficiencyProvider::RunA,(double)LUMINOSITYForRuns[i]/1000.);
        if(i==1) weight_provider[i]->setLumi(TopTriggerEfficiencyProvider::RunB,(double)LUMINOSITYForRuns[i]/1000.);
        if(i==2) weight_provider[i]->setLumi(TopTriggerEfficiencyProvider::RunC,(double)LUMINOSITYForRuns[i]/1000.);
        if(i==3) weight_provider[i]->setLumi(TopTriggerEfficiencyProvider::RunD,(double)LUMINOSITYForRuns[i]/1000.);
        if(i==4) {
            weight_provider[i]->setLumi(TopTriggerEfficiencyProvider::RunA,(double)LUMINOSITYForRuns[0]/1000.);
            weight_provider[i]->setLumi(TopTriggerEfficiencyProvider::RunB,(double)LUMINOSITYForRuns[1]/1000.);
            weight_provider[i]->setLumi(TopTriggerEfficiencyProvider::RunC,(double)LUMINOSITYForRuns[2]/1000.);
            weight_provider[i]->setLumi(TopTriggerEfficiencyProvider::RunD,(double)LUMINOSITYForRuns[3]/1000.);
        }
    }

    string PU_status="";
	switch(RunStatus_){
		case UncPUPlus:
            PU_status="_up";
			break;
		case UncPUMinus:
            PU_status="_down";
			break;
		default :
            PU_status="";
			break;
	}

    for(int i=0;i<NRuns;i++){
            sprintf(buffer_1,"PileUpReweight%s",RunsName[i].c_str());
            PileUpReweight[i] = new TH1D(buffer_1,"",60,0,60);
            sprintf(buffer_1,"interface/MyDataPileupHistogram_Run%s_json%s.root",RunsName[i].c_str(),PU_status.c_str());
            if(i==4)
                sprintf(buffer_1,"interface/MyDataPileupHistogram_AllRuns_json%s.root",PU_status.c_str());
            PUReweighting(PileUpReweight[i],buffer_1);
    }


	gROOT->ProcessLine(".L interface/setTDRStyle.C");
	gROOT->ProcessLine("setTDRStyle()");
	gStyle->SetErrorX(0.5);

	string jsonfile = "interface/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt";
	MakeJsonMap(jsonfile);

	double weight[samples_order_size];    
	for(int i=0;i<samples_order_size;i++) {
		switch(RunStatus_){
			case UncXsecPlus:
				weight[i] = (1+SAMPLE[i].unc)*(SAMPLE[i].xsec*LUMINOSITY)/(SAMPLE[i].ngen);
				break;
			case UncXsecMinus:
				weight[i] = (1-SAMPLE[i].unc)*(SAMPLE[i].xsec*LUMINOSITY)/(SAMPLE[i].ngen);
				break;
			default:
				weight[i] = (SAMPLE[i].xsec*LUMINOSITY)/(SAMPLE[i].ngen);
				break;
		}
	}
	int MCTruthChannels_[samples_order_size][4]; // [][hardronic, semileptonic, dilepton, total]    only for ttbar or t*    
    for(int isa_ = 0;isa_ < samples_order_size;isa_++) 
        for(int ich_ = 0;ich_ < 4;ich_++) MCTruthChannels_[isa_][ich_] = 0;

	// histograms
	const int Nmerge = 5;
	TH1F *ExcitedQuarkPlots[samples_order_size][_plots_order_size];
	TH1F *ExcitedQuarkPlotsSwap[_plots_order_size][samples_order_size];
	TH1F *ExcitedQuarkPlotsMerge[_plots_order_size][Nmerge];
	TH2D *ExcitedQuark2DPlots[samples_order_size][_plots2D_order_size];
	TH2D *ExcitedQuark2DPlotsSwap[_plots2D_order_size][samples_order_size];
	TH2D *ExcitedQuark2DPlotsMerge[_plots2D_order_size][Nmerge];

	TCanvas *c[_plots_order_size+ _plots2D_order_size]; 
	TLegend *legend_nm[_plots_order_size+ _plots2D_order_size];
	THStack *Stacks1[_plots_order_size+ _plots2D_order_size];
	THStack *Stacks1_[_plots_order_size+ _plots2D_order_size];
	for(int i =0;i<_plots_order_size+ _plots2D_order_size;i++) {
		sprintf(buffer_1,"c%i",i);
		c[i]= new TCanvas(buffer_1,"",640,640);
		sprintf(buffer_1,"Stacks1_%i",i);
		Stacks1[i] = new THStack(buffer_1,"");
		sprintf(buffer_1,"Stacks1_%i_",i);
		Stacks1_[i] = new THStack(buffer_1,"");
	}

    TFile *file_ = new TFile("MCTemplatesTree.root","recreate");

	for(int i=0;i<samples_order_size;i++) {

		std::cout<<"file : "<<SAMPLE[i].tag<<std::endl;
		for(int j=0;j<_plots_order_size;j++) {
			sprintf(buffer_1,"ExcitedQuarkPlots_%i_%i",i,j);
			ExcitedQuarkPlots[i][j] = new TH1F(buffer_1,"",PLOTS[j].Ndivide_,PLOTS[j].XminBond_,PLOTS[j].XmaxBond_);
		}

		for(int j=0;j<_plots2D_order_size;j++) {
			sprintf(buffer_1,"ExcitedQuark2DPlots_%i_%i",i,j);
			ExcitedQuark2DPlots[i][j] = new TH2D(buffer_1,"",
					PLOTS2D[j].XNdivide_,PLOTS2D[j].XminBond_,PLOTS2D[j].XmaxBond_,
					PLOTS2D[j].YNdivide_,PLOTS2D[j].YminBond_,PLOTS2D[j].YmaxBond_);
		}

        if( RunStatus_!= UncQsquare && 
                (i>=TTJets_matchingdown_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1 
                 && i<=DYJetsToLL_M_50_scaleup_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1))
            continue;
        if( !RunForMoreStat && 
                (i>=TTJets_FullLeptMGDecays_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v2&&
                i<=ZH_HToZZTo4L_M_126_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7C_v1))
            continue;

		analyze(i, ExcitedQuarkPlots[i], ExcitedQuark2DPlots[i] , weight[i], MCTruthChannels_[i]);

	}
    // Store information to get Efficiency(ISO) for weighting M(t+pho) without ISO
    for(int i=0;i<Nratio;i++){// [photon, ISR/FSR, gluon, light quark, heavy quark, electron, muon]
        hMCBGpTandEtaNumerator[i]->Write();
        hMCBGpTandEtaDenominator[i]->Write();
        hMCBGpTNumerator[i]->Write();
        hMCBGpTDenominator[i]->Write();
        hMCBGEtaNumerator[i]->Write();
        hMCBGEtaDenominator[i]->Write();
    }
    // Save for 1D template
	for(int i=0;i<_plots_order_size;i++) {
		for(int j=0;j<samples_order_size;j++) {
			sprintf(buffer_1,"templates/h%s%s.root",SAMPLE[j].tag,PLOTS[i].PlotName);
            ExcitedQuarkPlots[j][i]->SaveAs(buffer_1);
        }
    }
    //file_->Close();

	for(int i=0;i<_plots_order_size;i++) {
		for(int j=0;j<samples_order_size;j++) {
			ExcitedQuarkPlotsSwap[i][j] = (TH1F*) ExcitedQuarkPlots[j][i]->Clone();
		}
	}

	for(int i=0;i<_plots_order_size;i++) {
		for(int j=0;j<Nmerge;j++) {
			sprintf(buffer_1,"ExcitedQuarkPlotsMerge_%i_%i",i,j);
			ExcitedQuarkPlotsMerge[i][j] = new TH1F(buffer_1,"",PLOTS[i].Ndivide_,PLOTS[i].XminBond_,PLOTS[i].XmaxBond_);
		}
	}

    // Get Main Sysmatics uncertainty ( XSEC + QSCALE + MATCHING)
	for(int i=0;i<_plots_order_size;i++) {
        // false for loose and signal region, but true for control region (due to limited stats.)
        IncludeDYMatchScaleUnc = false;
        if(
               ( (i>=_N1p1lep5jet&&i<=_N2lep5jetRunD)||(i>=_MtstarChi2_1r1l5j&&i<=_MtstarChi2_0r2l5j)) &&
               (i!=_MtstarChi2_0r2l5j)

               )
            IncludeDYMatchScaleUnc = true;

        hMainSysUnc_[i] = MainSysUnc(ExcitedQuarkPlotsSwap[i], PLOTS[i].rebin_);
        sprintf(buffer_1,"hMainSysUnc_%i",i);
        hMainSysUnc_[i]->SetTitle(buffer_1);
    }
    /*
	for(int i=0;i<_plots2D_order_size;i++) {
		for(int j=0;j<samples_order_size;j++) {
			ExcitedQuark2DPlotsSwap[i][j] = (TH2D*) ExcitedQuark2DPlots[j][i]->Clone();
			std::cout<<PLOTS2D[i].PlotName<<" "<<SAMPLE[j].tag<<" : "<<ExcitedQuark2DPlotsSwap[i][j]->Integral() <<" ( "<<weight[j]<<" )"<<std::endl;
		}
	}
    */

	for(int i=0;i<_plots2D_order_size;i++) {
		for(int j=0;j<Nmerge;j++) {
			sprintf(buffer_1,"ExcitedQuark2DPlotsMerge_%i_%i",i,j);
			ExcitedQuark2DPlotsMerge[i][j] = new TH2D(buffer_1,"",
					PLOTS2D[i].XNdivide_,PLOTS2D[i].XminBond_,PLOTS2D[i].XmaxBond_,
					PLOTS2D[i].YNdivide_,PLOTS2D[i].YminBond_,PLOTS2D[i].YmaxBond_
					);
		}
	}
	// 2D template
	for(int i=0;i<_plots2D_order_size;i++){
		int template_indx = i;    // default = _ExcitedQuarkMass, test : _MR , _leading1jet, _leading2jets, _MfromMR
		sprintf(buffer_1,"hData%s",PLOTS2D[template_indx].PlotName);
		TH2D *hData = new TH2D(buffer_1,"",
				PLOTS2D[i].XNdivide_,PLOTS2D[i].XminBond_,PLOTS2D[i].XmaxBond_,
				PLOTS2D[i].YNdivide_,PLOTS2D[i].YminBond_,PLOTS2D[i].YmaxBond_);

		int Nsignal_=TGluonTgamma1500GeV-Tgamma500GeV+1;
		TH2D *hSigMCs[Nsignal_];
		for(int h_=0;h_<Nsignal_;h_++){
			sprintf(buffer_1,"hSigMC_%i%s",450+50*h_,PLOTS2D[template_indx].PlotName);
			hSigMCs[h_] = new TH2D(buffer_1,"",
					PLOTS2D[i].XNdivide_,PLOTS2D[i].XminBond_,PLOTS2D[i].XmaxBond_,
					PLOTS2D[i].YNdivide_,PLOTS2D[i].YminBond_,PLOTS2D[i].YmaxBond_);
			hSigMCs[h_]->Add(ExcitedQuark2DPlots[Tgamma500GeV+h_][template_indx]);
		}
		sprintf(buffer_1,"hMCBG%s",PLOTS2D[template_indx].PlotName);
		TH2D *hMCBG = new TH2D(buffer_1,"",
				PLOTS2D[i].XNdivide_,PLOTS2D[i].XminBond_,PLOTS2D[i].XmaxBond_,
				PLOTS2D[i].YNdivide_,PLOTS2D[i].YminBond_,PLOTS2D[i].YmaxBond_);
		hData->Add(ExcitedQuark2DPlots[Data][template_indx]);
		for(int i= TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
				i<=TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1;i++) 
			hMCBG->Add(ExcitedQuark2DPlots[i][template_indx]);
		hData->Rebin2D(PLOTS2D[template_indx].Xrebin_,PLOTS2D[template_indx].Yrebin_);
		for(int h_=0;h_<Nsignal_;h_++)
			hSigMCs[h_]->Rebin2D(PLOTS2D[template_indx].Xrebin_,PLOTS2D[template_indx].Yrebin_);
		hMCBG->Rebin2D(PLOTS2D[template_indx].Xrebin_,PLOTS2D[template_indx].Yrebin_);
		sprintf(buffer_1,"templates/hData%s.root",PLOTS2D[template_indx].PlotName);
		hData->SaveAs(buffer_1);
		for(int h_=0;h_<Nsignal_;h_++){
			sprintf(buffer_1,"templates/hSigMC_%i%s.root",450+50*h_,PLOTS2D[template_indx].PlotName);
			hSigMCs[h_]->SaveAs(buffer_1);
		}
		sprintf(buffer_1,"templates/hMCBG%s.root",PLOTS2D[template_indx].PlotName);
		hMCBG->SaveAs(buffer_1);
	}

	// For templates    -- start --
	// templates for limit
	//_MR2p2j,
	//_MtstarChi2,
	int template_indx = _N1p1lep4jet;    // default = _ExcitedQuarkMass, test : _MR , _leading1jet, _leading2jets, _MfromMR
	sprintf(buffer_1,"hData%s",RunStatusNames[RunStatus_].c_str());
	TH1F *hData = new TH1F(buffer_1,"",
			PLOTS[template_indx].Ndivide_,PLOTS[template_indx].XminBond_,PLOTS[template_indx].XmaxBond_); 


	int Nsignal_=TGluonTgamma1500GeV-Tgamma500GeV+1;
	TH1F *hSigMCs[Nsignal_];
	for(int h_=0;h_<Nsignal_;h_++){
		sprintf(buffer_1,"hSigMC_%i%s",450+50*h_,RunStatusNames[RunStatus_].c_str());
		hSigMCs[h_] = new TH1F(buffer_1,"",
				PLOTS[template_indx].Ndivide_,PLOTS[template_indx].XminBond_,PLOTS[template_indx].XmaxBond_); 
		hSigMCs[h_]->Add(ExcitedQuarkPlots[Tgamma500GeV+h_][template_indx]);
	}

	sprintf(buffer_1,"hMCBG%s",RunStatusNames[RunStatus_].c_str());
	TH1F *hMCBG = new TH1F(buffer_1,"",
			PLOTS[template_indx].Ndivide_,PLOTS[template_indx].XminBond_,PLOTS[template_indx].XmaxBond_); 
	hData->Add(ExcitedQuarkPlots[Data][template_indx]);
	//for(int i= QCD_Pt_15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2;
	//        i<=TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1;i++) 
	for(int i= TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
			i<=TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1;i++) 
		hMCBG->Add(ExcitedQuarkPlots[i][template_indx]);
	hData->Rebin(PLOTS[template_indx].rebin_);
	for(int h_=0;h_<Nsignal_;h_++)
		hSigMCs[h_]->Rebin(PLOTS[template_indx].rebin_);
	hMCBG->Rebin(PLOTS[template_indx].rebin_);

    TTree *ResultTree = new TTree("ResultTree","ResultTree");
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
    int ResultTree_RunMode = RunStatus_;
    ResultTree->Branch("ResultTree_sampleSize",&ResultTree_sampleSize,"ResultTree_sampleSize/I");
    ResultTree->Branch("ResultTree_plot",&ResultTree_plot,"ResultTree_plot/I");
    ResultTree->Branch("ResultTree_sampleID",&ResultTree_sampleID[0],"ResultTree_sampleID[ResultTree_sampleSize]/I");
    ResultTree->Branch("ResultTree_yield",&ResultTree_yield[0],"ResultTree_yield[ResultTree_sampleSize]/D");
    ResultTree->Branch("ResultTree_stat_error",&ResultTree_stat_error[0],"ResultTree_stat_error[ResultTree_sampleSize]/D");
    ResultTree->Branch("ResultTree_RunMode",&ResultTree_RunMode,"ResultTree_RunMode/I");

	for(int j=0;j<_plots_order_size;j++){
		std::cout<<"-------------"<<PLOTS[j].PlotName <<"-------------"<<std::endl;
		double NMC_=0;
		for(int i=QCD_Pt_15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2;i<samples_order_size;i++) 
			NMC_+=ExcitedQuarkPlots[i][j]->Integral();
		double NMCErr_=0;
		for(int i=QCD_Pt_15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2;i<samples_order_size;i++) {
			double err_=0;
			ExcitedQuarkPlots[i][j]->IntegralAndError(1,ExcitedQuarkPlots[i][j]->GetNbinsX(),err_);
			NMCErr_+=err_*err_;
		}
		std::cout<<"Total MC : "<<NMC_<<" +/- "<< sqrt(NMCErr_)<<std::endl;
        ResultTree_plot = j;
        ResultTree_sampleID[0] = -1;
        ResultTree_yield[0] = NMC_;
        ResultTree_stat_error[0] = sqrt(NMCErr_);
		for(int i=0;i<samples_order_size;i++){ 
			double err_=0.;
			double yield_ = ExcitedQuarkPlots[i][j]->IntegralAndError(1,ExcitedQuarkPlots[i][j]->GetNbinsX(),err_);
			if(err_==0) err_ = weight[i];
			std::cout<<SAMPLE[i].tag<<" : "<<
				yield_ <<
				" +/- "<< err_<<std::endl;
			//std::cout<<SAMPLE[i].tag<<" : "<<ExcitedQuarkPlots[i][j]->Integral() <<" ( "<<weight[i]<<" )"<<std::endl;

            ResultTree_sampleID[i+1] = i;
            ResultTree_yield[i+1] = yield_;
            ResultTree_stat_error[i+1] = err_;
		}
        ResultTree->Fill();
	}
    ResultTree->Write();

    int sampleForMCTruthChannel = TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1;
    std::cout<<SAMPLE[sampleForMCTruthChannel].tag<<" (had, semi, dilep, total) = ("<<
        MCTruthChannels_[sampleForMCTruthChannel][0]<<" , "<<
        MCTruthChannels_[sampleForMCTruthChannel][1] << " , "<<
        MCTruthChannels_[sampleForMCTruthChannel][2] << " , "<<
        MCTruthChannels_[sampleForMCTruthChannel][3]<<" ) : "
            << 100.* MCTruthChannels_[sampleForMCTruthChannel][1]/ 
               (MCTruthChannels_[sampleForMCTruthChannel][0]+ 
                MCTruthChannels_[sampleForMCTruthChannel][1]+
                MCTruthChannels_[sampleForMCTruthChannel][2])<<" %"<<std::endl;

    // [e_s, e_l, e_j, eL_l, eL_j]
    string RatioNames[5] = {
        "e_s",
        "e_l",
        "e_j",
        "eL_l",
        "eL_j"
    };


    int _PdgID = _N2p1lep4jet;
    /*
	for(int i=0;i<samples_order_size;i++) {
		sprintf(buffer_1,"templates/hPdgID_%s.root",SAMPLE[i].tag);
		ExcitedQuarkPlots[i][_PdgID]->SaveAs(buffer_1);
		// all + 18 (==0)
		//   (r,g/q,l,z,w,unknown) = 
		//   (22,21/1-6 or -1~-6,11~15 or -11~-15, 23, 24, 0)

		std::cout<<SAMPLE[i].tag<<"(r,g/q,l,z,w,unknown) : "<<
			ExcitedQuarkPlots[i][_PdgID]->Integral(22+18,22+18)/ExcitedQuarkPlots[i][_PdgID]->Integral() *100 <<" "<<
			(ExcitedQuarkPlots[i][_PdgID]->Integral(-6+18,6+18) - ExcitedQuarkPlots[i][_PdgID]->Integral(18,18) + ExcitedQuarkPlots[i][_PdgID]->Integral(21+18,21+18) )/ExcitedQuarkPlots[i][_PdgID]->Integral() *100 <<" "<<
(ExcitedQuarkPlots[i][_PdgID]->Integral(11+18,15+18)+ExcitedQuarkPlots[i][_PdgID]->Integral(-15+18,-11+18))/ExcitedQuarkPlots[i][_PdgID]->Integral() *100 <<" "<<
	ExcitedQuarkPlots[i][_PdgID]->Integral(23+18,23+18)/ExcitedQuarkPlots[i][_PdgID]->Integral() *100 <<" "<<
	ExcitedQuarkPlots[i][_PdgID]->Integral(24+18,24+18)/ExcitedQuarkPlots[i][_PdgID]->Integral() *100 <<" "<<
	ExcitedQuarkPlots[i][_PdgID]->Integral(18,18)/ExcitedQuarkPlots[i][_PdgID]->Integral() *100 <<" "
													 <<std::endl;
	}
    */
	std::cout<<"---------------------------------------"<<std::endl;

	hData->SaveAs("templates/hData.root");

	for(int h_=0;h_<Nsignal_;h_++){
		sprintf(buffer_1,"templates/hSigMC_%i%s.root",450+50*h_,RunStatusNames[RunStatus_].c_str());
		hSigMCs[h_]->SaveAs(buffer_1);
	}

	sprintf(buffer_1,"templates/hMCBG%s.root",RunStatusNames[RunStatus_].c_str());
	hMCBG->SaveAs(buffer_1);
	// For templates    -- end --

    string plotsName_ = "";
	for(int i =0;i<_plots_order_size+_plots2D_order_size;i++) {
		c[i]->cd();

		if(i<_plots_order_size){
			HistMerge(ExcitedQuarkPlotsSwap[i], ExcitedQuarkPlotsMerge[i]);
			DRAWSTACK(c[i],Stacks1[i],
					ExcitedQuarkPlotsMerge[i],
					ExcitedQuarkPlots[Tgamma950GeV][i],
					ExcitedQuarkPlots[Data][i],
					legend_nm[i],
					PLOTS[i].YaxisTitle_,PLOTS[i].XaxisTitle_,PLOTS[i].XUnit_, PLOTS[i].rebin_, hMainSysUnc_[i]);
            plotsName_ = PLOTS[i].PlotName;
		}else{
            /*
			HistMerge(ExcitedQuark2DPlotsSwap[i-_plots_order_size], ExcitedQuark2DPlotsMerge[i-_plots_order_size]);
			DRAWSTACK(Stacks1[i],
					ExcitedQuark2DPlotsMerge[i-_plots_order_size],
					ExcitedQuark2DPlots[Tgamma950GeV][i-_plots_order_size],
					ExcitedQuark2DPlots[Data][i-_plots_order_size],
					legend_nm[i],
					PLOTS2D[i-_plots_order_size].YaxisTitle_,PLOTS2D[i-_plots_order_size].XaxisTitle_,
					PLOTS2D[i-_plots_order_size].XUnit_, PLOTS2D[i-_plots_order_size].Xrebin_,
					PLOTS2D[i-_plots_order_size].YUnit_, PLOTS2D[i-_plots_order_size].Yrebin_
				 );
            */
            plotsName_ = PLOTS2D[i-_plots_order_size].PlotName;
		}

        sprintf(buffer_1,"PLOTS_comv2/%s.pdf",plotsName_.c_str());
        c[i]->SaveAs(buffer_1);
        sprintf(buffer_1,"PLOTS_comv2/%s.png",plotsName_.c_str());
        c[i]->SaveAs(buffer_1);
        sprintf(buffer_1,"PLOTS_comv2/%s.root",plotsName_.c_str());
        c[i]->SaveAs(buffer_1);

	}


	TCanvas *checkpdf = new TCanvas("checkpdf","",640,640);
	checkpdf->cd();
	for(int i =0;i<_plots_order_size+_plots2D_order_size;i++) {
		if(i<_plots_order_size){
			DRAWSTACK(c[i],Stacks1_[i],
					ExcitedQuarkPlotsMerge[i],
					ExcitedQuarkPlots[Tgamma950GeV][i],
					ExcitedQuarkPlots[Data][i],
					legend_nm[i],
					PLOTS[i].YaxisTitle_,PLOTS[i].XaxisTitle_,PLOTS[i].XUnit_, 1,hMainSysUnc_[i]);
		}else{
            /*
			DRAWSTACK(Stacks1[i],
					ExcitedQuark2DPlotsMerge[i-_plots_order_size],
					ExcitedQuark2DPlots[Tgamma950GeV][i-_plots_order_size],
					ExcitedQuark2DPlots[Data][i-_plots_order_size],
					legend_nm[i],
					PLOTS2D[i-_plots_order_size].YaxisTitle_,PLOTS2D[i-_plots_order_size].XaxisTitle_,
					PLOTS2D[i-_plots_order_size].XUnit_, PLOTS2D[i-_plots_order_size].Xrebin_,
					PLOTS2D[i-_plots_order_size].YUnit_, PLOTS2D[i-_plots_order_size].Yrebin_
				 );
            */
		}
		if(_plots_order_size+_plots2D_order_size!=1){
			if(i==0){
				checkpdf->Print("checkpdf_comv2.pdf(");
			}else if(i==_plots_order_size+_plots2D_order_size-1){
				checkpdf->Print("checkpdf_comv2.pdf)");
			}else{
				checkpdf->Print("checkpdf_comv2.pdf");
			}
		}else{
			checkpdf->Print("checkpdf_comv2.pdf");
		}
	}

	delete checkpdf;

	time2=time(NULL); time3=time2-time1; cout<<"---- Spending time : " << (time3/60.0) <<" (mins) ----"<<endl;
	std::cout<<"done "<<std::endl;
}
