#ifndef ModelParameter_H
#define ModelParameter_H
#include <cmath>

#define Pi 3.1415926535897932384626433832795029
#define Pi2 (Pi*Pi)
#define PiHalf (Pi/2)

#define MZ 91.1876
#define MZ2 (MZ*MZ)
#define MW 80.385
#define MW2 (MW*MW)
#define CW (MW/MZ)
#define w (acos(CW))
#define SW (sin(w))
#define SW2 (SW*SW)
#define CW2 (CW*CW)

#define GF 1.16637e-5
#define Alfa0 (1/137.035999074)
#define AlfaGF (sqrt2/pi*GF*MW2*SW2)
#define AlfaMZ (1/127.944)
#define Alfa Alfa0
#define EL (sqrt(4*Pi*Alfa))
#define Alfa2 (Alfa*Alfa)
#define DeltaAlfa5Had .027547
#define AlfasMZ .1184
#define Alfas 0.0
#define Alfashgg 0.1184

#define vev (2.0*MW*SW/EL)
#define vev2 (vev*vev)


#define Mh1 125.09
#define Mh12 (Mh1*Mh1)

#define ME .5109989280e-3
#define ME2 (ME*ME)
#define MM 105.6583715e-3
#define MM2 (MM*MM)
#define ML 1776.82e-3
#define ML2 (ML*ML)


#define MU 7.356e-2
#define MU2 (MU*MU)
#define MC 1.275
#define MC2 (MC*MC)
#define MT 173.21
#define MT2 (MT*MT)

#define MD MU
#define MD2 (MD*MD)
#define MS 95e-3
#define MS2 (MS*MS)
#define MB1S 4.66
#define MBatMB 4.18
#define MB 4.66
#define MB2 (MB*MB)

#define Sin(i) sin(i)
#define Csc(i) (1/sin(i))
#define Cos(i) cos(i)
#define Sec(i) (1/cos(i))
#define Tan(i) tan(i)
#define Cot(i) (1/tan(i))
#define Power(i,j) (pow(i,j))


bool CheckUnitarity(double Mh22, double MHA2, double theta, double vs, double d2);
bool CheckStability(double Mh22, double MHA2, double theta, double vs, double d2);

#endif