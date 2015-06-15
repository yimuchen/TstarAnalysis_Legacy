#ifndef MCTruthChannel_H
#define MCTruthChannel_H

#include "format.h"
#include <iostream>
/*
// Only for ttbar
Wlight -> ud, us
Wheavy -> jj (excluding ud, us)
Wleptonic -> lv

top topbar decay channels :

Full hadronic decay
0). Wlight + Wheavy
1). Wlight + Wlight
2). Wheavy + Wheavy

Semi-leptonic decay
3). Wlight + Wleptonic
4). Wheavy + Wleptonic

Double-leptonic decay
5). Wleptonic + Wleptonic
 */

int MCTruthChannel(GenInfoBranches GenInfo){
    int indexWp_=-1;
    int indexWm_=-1;
    for(int i=0;i<GenInfo.Size;i++){
        if(GenInfo.PdgID[i]==24){
            indexWp_ = i;
        }else if(GenInfo.PdgID[i]==-24){
            indexWm_ = i;
        }
    }
    if (indexWp_==-1|| indexWm_==-1) return -1;

    bool Wplight = fabs(GenInfo.PdgID[GenInfo.Da1[indexWp_]])<=3 && fabs(GenInfo.PdgID[GenInfo.Da2[indexWp_]])<=3;
    bool Wpheavy = !(fabs(GenInfo.PdgID[GenInfo.Da1[indexWp_]])<=3 && fabs(GenInfo.PdgID[GenInfo.Da2[indexWp_]])<=3) && 
        fabs(GenInfo.PdgID[GenInfo.Da1[indexWp_]])< 10;
    bool Wplepton= fabs(GenInfo.PdgID[GenInfo.Da1[indexWp_]])>10;

    bool Wmlight = fabs(GenInfo.PdgID[GenInfo.Da1[indexWm_]])<=3 && fabs(GenInfo.PdgID[GenInfo.Da2[indexWm_]])<=3;
    bool Wmheavy = !(fabs(GenInfo.PdgID[GenInfo.Da1[indexWm_]])<=3 && fabs(GenInfo.PdgID[GenInfo.Da2[indexWm_]])<=3) && 
        fabs(GenInfo.PdgID[GenInfo.Da1[indexWm_]])< 10;
    bool Wmlepton= fabs(GenInfo.PdgID[GenInfo.Da1[indexWm_]])>10;

    if((Wplight&&Wmheavy)||(Wpheavy&&Wmlight)){
        return 0;
    }else if(Wplight&&Wmlight){
        return 0;
        //return 1;
    }else if(Wpheavy&&Wmheavy){
        return 0;
        //return 2;
    }else if((Wplight&&Wmlepton)||(Wplepton&&Wmlight)){
        //return 3;
        return 1;
    }else if((Wpheavy&&Wmlepton)||(Wplepton&&Wmheavy)){
        //return 4;
        return 1;
    }else if(Wplepton && Wmlepton){
        //return 5;
        return 2;
    }else{
        std::cout<<"[WARNING] : This event is not considered in the MC truth"<<std::endl;
    }
}

#endif
