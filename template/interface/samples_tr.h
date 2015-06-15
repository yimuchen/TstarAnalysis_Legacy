#ifndef SAMPLES_tr_H
#define SAMPLES_tr_H

#include "ConstantNumbers.h"

enum samples_tr_order{
    Data,
    Tgamma450GeV,
    Tgamma500GeV,
    Tgamma550GeV,
    Tgamma600GeV,
    Tgamma650GeV,
    Tgamma700GeV,
    Tgamma750GeV,
    Tgamma800GeV,
    Tgamma850GeV,
    Tgamma950GeV,
    Tgamma1100GeV,
    Tgamma1300GeV,
    QCD_Pt_15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2,
    QCD_Pt_30to50_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2,
    QCD_Pt_170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_50to80_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    QCD_Pt_80to120_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v2,
    QCD_Pt_800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTJets_MassiveBinDECAY_TuneZ2star_8TeV_madgraph_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    T_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    T_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v2,
    T_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    Tbar_s_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    Tbar_t_channel_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    Tbar_tW_channel_DR_TuneZ2star_8TeV_powheg_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DYJetsToLL_M_50_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v1,
    DiPhotonJets_8TeV_madgraph_tarball_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WJetsToLNu_TuneZ2Star_8TeV_madgraph_tarball_Summer12_DR53X_PU_S10_START53_V7A_v2,
    WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X_PU_S10_START53_V7A_v1,
    //samesignWWjj_Summer11MG7TeV_AlexisLHE,
    TTWJets_8TeV_madgraph_Summer12_DR53X_PU_S10_START53_V7A_v1,
    TTZJets_8TeV_madgraph_v2_Summer12_DR53X_PU_S10_START53_V7A_v1,
    samples_order_size
};

struct sample_t{
    char filename[256];
    char tag[64];
    double ngen;
    double xsec;
    double unc; // 1 : 100%, 0.11 : 11%, 0 : no uncertainty
};


struct sample_t SAMPLE[samples_order_size] = {
    {"./REDUCE_DATA/PhotonData/diphotondijet_*_Run2012*.root",
        "data", LUMINOSITY , 1 , 0},
    {"./REDUCE_DATA/diphotondijet_Tgamma450GeV.root",
        "Tgamma450GeV", 28546 , 11.984291 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/diphotondijet_Tgamma500GeV.root",
        "Tgamma500GeV", 29244 , 5.087496 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/diphotondijet_Tgamma550GeV.root",
        "Tgamma550GeV", 18898 , 2.315665 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/diphotondijet_Tgamma600GeV.root",
        "Tgamma600GeV", 20397 , 1.115097 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/diphotondijet_Tgamma650GeV.root",
        "Tgamma650GeV", 12728 , 0.562451 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/diphotondijet_Tgamma700GeV.root",
        "Tgamma700GeV", 16898 , 0.294853 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/diphotondijet_Tgamma750GeV.root",
        "Tgamma750GeV", 10883 , 0.159653 , 0},  // NNLO from excited quark
    {"./REDUCE_DATA/diphotondijet_Tgamma800GeV.root",
        "Tgamma800GeV", 38966 , 0.088853 , 0},  // NNLO from excited quark
        //"Tgamma800GeV", 38966 , 0.0208 , 0},  // NNLO from HATHOR
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-15to30_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "QCD_Pt-15to30", 361461 , 988287420 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-30to50_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-30to50", 413046 , 66285328 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-120to170_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "QCD_Pt-120to170", 883691 , 156293.3 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-170to300_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-170to300", 295728 , 34138.15 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-300to470_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-300to470", 315004 , 1759.549 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-470to600_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-470to600", 324593 , 113.8791 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-50to80_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-50to80", 385138 , 8148778.0 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-600to800_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-600to800", 265508 , 26.9921 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-80to120_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "QCD_Pt-80to120", 844107 , 1033680.0 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "QCD_Pt-800to1000", 313083 , 3.550036 , 1},  // LO
    {"./REDUCE_DATA/diphotondijet_TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "TTJets", 6759380 , 227 , 0.11},    // TOP-12-007
    {"./REDUCE_DATA/diphotondijet_T_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "T_TuneZ2_s-channel", 259575 , 3.79 , 0.3}, // NNLO
    {"./REDUCE_DATA/diphotondijet_T_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "T_TuneZ2_t-channel", 23745 , 56.4 , 0.3},  // NNLO
    {"./REDUCE_DATA/diphotondijet_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "T_TuneZ2_tW-channel", 496681 , 11.1 , 0.3},    // NNLO
    {"./REDUCE_DATA/diphotondijet_Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tbar_TuneZ2_s-channel", 139803 , 1.76 , 0.3},  // NNLO
    {"./REDUCE_DATA/diphotondijet_Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tbar_TuneZ2_t-channel", 1932776 , 30.7 , 0.3}, // NNLO
    {"./REDUCE_DATA/diphotondijet_Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "Tbar_TuneZ2_tW-channel", 492545 , 11.1 , 0.3}, // NNLO
    {"./REDUCE_DATA/diphotondijet_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "DYJetsToLL_TuneZ2_M-50", 22491289  , 1120.*3*3503.71/3351.97 , 0.04},  // SMP-12-011
    {"./REDUCE_DATA/diphotondijet_DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "DiPhotonJets", 1152621 , 75.39 ,0.},   // CTEQ
    {"./REDUCE_DATA/diphotondijet_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v2.root",
        "WJetsToLNu", 54643853 , 36257.2 , 0.03},   // NNLO
    {"./REDUCE_DATA/diphotondijet_WW_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "WW", 9929882 , 54.84 , 0.26},  // CTEQ
    {"./REDUCE_DATA/diphotondijet_WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "WZ", 9989440 , 33.21 ,0.30},  // CTEQ
    {"./REDUCE_DATA/diphotondijet_ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ZZ", 9789109 , 8.059 ,0.21},   // CTEQ
    //{"./REDUCE_DATA/diphotondijet_samesignWWjj_Summer11MG7TeV_AlexisLHE.root",
    //    "samesignWWjj", 51600 , 0.159 ,0.49},
    {"./REDUCE_DATA/diphotondijet_TTWJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ttW", 105604 , 0.2149 , 0.19},
    {"./REDUCE_DATA/diphotondijet_TTZJets_8TeV-madgraph_v2_Summer12_DR53X-PU_S10_START53_V7A-v1.root",
        "ttZ", 209385 , 0.172 ,0.28}
};

#endif
