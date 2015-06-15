#ifndef ConstantNumbers_H
#define ConstantNumbers_H

//const float LUMINOSITY       = 1421.2; 
//const float LUMINOSITY       = 16600; 
const float LUMINOSITY       = 19625; 
const double ELECTRON_MASS   = 0.0005109989;
const double MUON_MASS       = 0.105658;
const double Z_MASS          = 91.1876;
const double W_MASS          = 80.398;
const double b_MASS          = 4.19;
const double top_MASS        = 172.9;
const double dRlrOffset      = 0.01;

enum RunStatus{
    Normal,
    UncXsecPlus,
    UncXsecMinus,
    UncPUPlus,
    UncPUMinus,
    RunStatusSize
};

string RunStatusNames[RunStatusSize] = {
    "",
    "UncXsecPlus",
    "UncXsecMinus",
    "UncPUPlus",
    "UncPUMinus"
};

#endif
