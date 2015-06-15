#ifndef EffectiveAreaElectron_H
#define EffectiveAreaElectron_H



float EffectiveAreaElectron(float eta){
	// https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection
	float EA[7] = {
        0.13,
        0.14,
        0.07,
        0.09,
        0.11,
        0.11,
        0.14
	};



	if (fabs(eta)<1.0)
		return EA[0];
	if (1.0<=fabs(eta)&&fabs(eta)<1.479)
		return EA[1];
	if (1.479<=fabs(eta)&&fabs(eta)<2.0)
		return EA[2];
	if (2.0<=fabs(eta)&&fabs(eta)<2.2)
		return EA[3];
	if (2.2<=fabs(eta)&&fabs(eta)<2.3)
		return EA[4];
	if (2.3<=fabs(eta)&&fabs(eta)<2.4)
		return EA[5];
	if (fabs(eta)>2.4)
		return EA[6];

}


#endif
