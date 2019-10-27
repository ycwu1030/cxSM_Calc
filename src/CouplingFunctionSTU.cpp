#include "CouplingFunctionSTU.h"
#include "ModelParameters.h"

namespace STU{
ComplexType SinSMCS(double Mh22, double MHA2, double theta, double vs, double d2)
{

	ComplexType FreeTerm = ((-(MZ2*B0i(bb0,0,Mh12,MZ2)) + MZ2*B0i(bb0,0,Mh22,MZ2) + MZ2*B0i(bb0,MZ2,Mh12,MZ2) - MZ2*B0i(bb0,MZ2,Mh22,MZ2) + B0i(bb00,0,Mh12,MZ2) - B0i(bb00,0,Mh22,MZ2) - B0i(bb00,MZ2,Mh12,MZ2) + B0i(bb00,MZ2,Mh22,MZ2))*Power(Sin(theta),2))/(MZ2*Pi);
	return FreeTerm;

}

ComplexType SinSMCS_Z2(double Mh22, double MHA2, double d2, double del2)
{

	ComplexType FreeTerm = 0;
	return FreeTerm;

}

ComplexType TinSMCS(double Mh22, double MHA2, double theta, double vs, double d2)
{

	ComplexType FreeTerm = ((MW2*B0i(bb0,0,Mh12,MW2) - MZ2*B0i(bb0,0,Mh12,MZ2) - MW2*B0i(bb0,0,Mh22,MW2) + MZ2*B0i(bb0,0,Mh22,MZ2) - B0i(bb00,0,Mh12,MW2) + B0i(bb00,0,Mh12,MZ2) + B0i(bb00,0,Mh22,MW2) - B0i(bb00,0,Mh22,MZ2))*Power(Sin(theta),2))/(4.*MW2*Pi*SW2);
	return FreeTerm;

}

ComplexType TinSMCS_Z2(double Mh22, double MHA2, double d2, double del2)
{

	ComplexType FreeTerm = 0;
	return FreeTerm;

}

ComplexType UinSMCS(double Mh22, double MHA2, double theta, double vs, double d2)
{

	ComplexType FreeTerm = -(((MW2*MZ2*B0i(bb0,0,Mh12,MW2) - MW2*MZ2*B0i(bb0,0,Mh12,MZ2) - MW2*MZ2*B0i(bb0,0,Mh22,MW2) + MW2*MZ2*B0i(bb0,0,Mh22,MZ2) - MW2*MZ2*B0i(bb0,MW2,Mh12,MW2) + MW2*MZ2*B0i(bb0,MW2,Mh22,MW2) + MW2*MZ2*B0i(bb0,MZ2,Mh12,MZ2) - MW2*MZ2*B0i(bb0,MZ2,Mh22,MZ2) - MZ2*B0i(bb00,0,Mh12,MW2) + MW2*B0i(bb00,0,Mh12,MZ2) + MZ2*B0i(bb00,0,Mh22,MW2) - MW2*B0i(bb00,0,Mh22,MZ2) + MZ2*B0i(bb00,MW2,Mh12,MW2) - MZ2*B0i(bb00,MW2,Mh22,MW2) - MW2*B0i(bb00,MZ2,Mh12,MZ2) + MW2*B0i(bb00,MZ2,Mh22,MZ2))*Power(Sin(theta),2))/(MW2*MZ2*Pi));
	return FreeTerm;

}

ComplexType UinSMCS_Z2(double Mh22, double MHA2, double d2, double del2)
{

	ComplexType FreeTerm = 0;
	return FreeTerm;

}

} //end namespace STU
