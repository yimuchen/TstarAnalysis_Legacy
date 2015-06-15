#ifndef ConstantNumbers_H
#define ConstantNumbers_H

#include <string>

using namespace std;
//const float LUMINOSITY       = 1421.2; 
//const float LUMINOSITY       = 16600; 
//const float LUMINOSITY       = 19625; 
const float LUMINOSITY       = 19712; 
const double ELECTRON_MASS   = 0.0005109989;
const double MUON_MASS       = 0.105658;
const double Z_MASS          = 91.1876;
const double W_MASS          = 80.385;//80.398;
const double b_MASS          = 4.18;//4.19;
const double top_MASS        = 173.21;//172.9;
const double dRlrOffset      = 0.01;

enum RunStatus{
    Normal,
    UncXsecPlus,
    UncXsecMinus,
    UncPUPlus,
    UncPUMinus,
    UncJESPlus,
    UncJESMinus,
    UncJERPlus,
    UncJERMinus,
    UncQsquare,
    UncTrigPlus,
    UncTrigMinus,
    UncTopPtPlus,
    UncTopPtMinus,
    UncPhoIDPlus,
    UncPhoIDMinus,
    UncLepIDPlus,
    UncLepIDMinus,
    RunStatusSize
};

string RunStatusNames[RunStatusSize] = {
    "",
    "UncXsecPlus",
    "UncXsecMinus",
    "UncPUPlus",
    "UncPUMinus",
    "UncJESPlus",
    "UncJESMinus",
    "UncJERPlus",
    "UncJERMinus",
    "UncQsquare",
    "UncTrigPlus",
    "UncTrigMinus",
    "UncTopPtPlus",
    "UncTopPtMinus",
    "UncPhoIDPlus",
    "UncPhoIDMinus",
    "UncLepIDPlus",
    "UncLepIDMinus"
};

enum Photon_WP {
    P_LOOSE,
    P_MEDIUM,
    P_TIGHT
};


#endif
