#include <iostream>
#include "samples_tr_1lepton4jets.h"


void xsec(){

    for(int i=0;i<samples_order_size;i++)
        printf("| %s | %f (pb) | %f (/pb) | \n",
                SAMPLE[i].tag,
                SAMPLE[i].xsec,
                SAMPLE[i].ngen/SAMPLE[i].xsec
                );

}
