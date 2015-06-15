#ifndef MCTRUTH_H
#define MCTRUTH_H 

#include "math.h"

/*
   enum Jet_Labels {
   isr_label = 0,
   lepb_label = 11,
   hadb_label = 12,
   hadw1_label = 13,
   hadw2_label = 14,
   higgs_label = 15,
   gluon1_label = 16,
   gluon2_label = 17,
   unknown_label = 20
   };

   For Excited Quark : 
   (lepb_label+Wlep+gluon1_label) = (hadb_label+hadw1_label+hadw2_label+gluon2_label)

   gluon2_label, gluon1_label, hadw2_label, hadw1_label, hadb_label, lepb_label
   2^5       ,  2^4      ,  2^3       ,  2^2       ,  2^1        ,  2^0
 */


int MCTruth(int JetType[],int JetGenJetType[], int N){

    int MCTruth_ = 0;

    for(int i=0;i<N;i++){
        if(JetType[i]==11&&(JetType[i]==JetGenJetType[i])) MCTruth_ += pow(2.,0);
        if(JetType[i]==12&&(JetType[i]==JetGenJetType[i])) MCTruth_ += pow(2.,1);
        if(JetType[i]==13&&(JetType[i]==JetGenJetType[i])) MCTruth_ += pow(2.,2);
        if(JetType[i]==14&&(JetType[i]==JetGenJetType[i])) MCTruth_ += pow(2.,3);
        if(JetType[i]==16&&(JetType[i]==JetGenJetType[i])) MCTruth_ += pow(2.,4);
        if(JetType[i]==17&&(JetType[i]==JetGenJetType[i])) MCTruth_ += pow(2.,5);
    }

    return MCTruth_;
}

#endif
