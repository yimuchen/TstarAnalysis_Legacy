#ifndef SolutionOfWNeutrino_H
#define SolutionOfWNeutrino_H

#include <iostream>
#include "TLorentzVector.h"
#include <math.h>

void SolutionOfWNeutrino(TLorentzVector vlep,float MET, float METPhi, float wmass, float& nuz1, float& nuz2)
//
// Purpose: Solve for the neutrino longitudinal z-momentum that makes
//          the leptonic W have mass WMASS.
//
// Inputs:
//   ev -          The event to solve.
//   wmass -       The desired W mass.
//
// Outputs:
//   nuz1 -        First solution (smaller absolute value).
//   nuz2 -        Second solution. 
//                           
// Returns:                  
//   True if there was a real solution.  False if there were only
//   imaginary solutions.  (In that case, we just set the imaginary
//   part to zero.)
{
    float vnu_x = MET*cos(METPhi);
    float vnu_y = MET*sin(METPhi);
    float vnu_perp2 = vnu_x*vnu_x + vnu_y*vnu_y;

    float x = vlep.Px()*vnu_x + vlep.Py()*vnu_y + wmass*wmass/2;
    float a = vlep.Pz()*vlep.Pz() - vlep.E()*vlep.E();
    float b = 2*x*vlep.Pz();
    float c = x*x - vnu_perp2 * vlep.E()*vlep.E();

    float d = b*b - 4*a*c;
    if (d < 0) {
        d = 0;
    }

    nuz1 = (-b + sqrt (d))/2/a;
    nuz2 = (-b - sqrt (d))/2/a;
    if (fabs (nuz1) > fabs (nuz2))
        std::swap (nuz1, nuz2);

}

#endif
