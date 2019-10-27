#ifndef CouplingFunctionSTU
#define CouplingFunctionSTU
#include "clooptools.h"

namespace STU{
ComplexType SinSMCS(double Mh22, double MHA2, double theta, double vs, double d2);

ComplexType SinSMCS_Z2(double Mh22, double MHA2, double d2, double del2);

ComplexType TinSMCS(double Mh22, double MHA2, double theta, double vs, double d2);

ComplexType TinSMCS_Z2(double Mh22, double MHA2, double d2, double del2);

ComplexType UinSMCS(double Mh22, double MHA2, double theta, double vs, double d2);

ComplexType UinSMCS_Z2(double Mh22, double MHA2, double d2, double del2);

} //end namespace STU
#endif
