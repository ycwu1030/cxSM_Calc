#include "CouplingFunctionSMCS.h"
#include "ModelParameters.h"

namespace SMCS{
ComplexType dKappahW(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(56 - 11*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(Power(EL,2)*(-2 + 23*Cos(2*w) + 30*Cos(4*w) - 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(13 + 11*Cos(2*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,m12,Mh12,Mh12))*(1.0*(-3*EL*Csc(w)*(MW*Power(Cos(theta),2)*(d2*Power(vs,2) + 2*(Mh12 - Mh22)*Power(Cos(theta),2))*Power(Sin(theta),3) + EL*Mh12*vs*Csc(w)*(-1 + Power(Cos(theta),7) + Power(Cos(theta),5)*Power(Sin(theta),2)))))/(1.0*(128*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm9 = (B0i(bb0,m12,Mh12,Mh22))*(1.0*(EL*Cos(theta)*Csc(w)*Sin(theta)*(-2*EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + 2*MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta)))/(1.0*(256*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm10 = (B0i(bb0,m12,Mh22,Mh22))*(1.0*(-(EL*Csc(w)*Power(Sin(theta),2)*(MW*(Mh12 - Mh22 + 6*d2*Power(vs,2))*Cos(theta) + 3*(Mh12 - Mh22)*MW*Cos(3*theta) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta))*Sin(2*theta))))/(1.0*(512*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm11 = (B0i(bb0,m12,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12 + 12*MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm12 = (B0i(bb0,m12,MZ2,MZ2))*(1.0*(Power(EL,2)*(Mh12 + 24*MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm13 = (B0i(bb0,m22,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm14 = (B0i(bb0,m22,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,m22,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,m22,MW2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*Power(Tan(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,m32,0,MW2))*(1.0*(5*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm18 = (B0i(bb0,m32,MB2,MT2))*(1.0*(-3*Power(EL,2)*(MB2 + MT2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(8*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm19 = (B0i(bb0,m32,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm20 = (B0i(bb0,m32,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm21 = (B0i(bb0,m32,MW2,MZ2))*(1.0*(-(Power(EL,2)*(23 + 36*Cos(2*w) + 5*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm24 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm25 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm26 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm27 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm29 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(-8*MW2*Power(Csc(2*w),2) + (Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(Mh12 - Mh22))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm31 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm32 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm33 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm34 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm35 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb0,MW2,0,MW2))*(1.0*(Power(EL,2)*(-1 + 5*Cos(2*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2*Cos(2*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(22 + 45*Cos(2*w) + 10*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(256*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),2))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(5 + 4*Cos(2*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm69 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm70 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm71 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm72 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm81 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm82 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm83 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm84 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm85 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm86 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm87 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm88 = (B0i(dbb0,MW2,0,MW2))*(1.0*(Power(EL,2)*MW2*(-1 + Cos(theta))))/(1.0*(4*Power(Pi,2)));
	ComplexType PVTerm89 = (B0i(dbb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-MT2 + MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm90 = (B0i(dbb0,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*MW2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm91 = (B0i(dbb0,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*MW2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm92 = (B0i(dbb0,MW2,MW2,MZ2))*(1.0*(-(Power(EL,2)*MW2*(7 + 12*Cos(2*w) + Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm93 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm94 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm95 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm96 = (B0i(dbb00,MW2,0,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(2*Power(Pi,2)));
	ComplexType PVTerm97 = (B0i(dbb00,MW2,MB2,MT2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm98 = (B0i(dbb00,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm99 = (B0i(dbb00,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm100 = (B0i(dbb00,MW2,MW2,MZ2))*(1.0*(-(Power(EL,2)*(5 + 4*Cos(2*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm101 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm102 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm103 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm104 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*Power(Sec(w),2))))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm105 = (B0i(dbb1,MW2,0,MW2))*(1.0*(Power(EL,2)*MW2*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm106 = (B0i(dbb1,MW2,MB2,MT2))*(1.0*(-3*Power(EL,2)*MW2*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm107 = (B0i(dbb1,MW2,MW2,MZ2))*(1.0*(-(Power(EL,2)*MW2*(-1 + Cos(theta))*Power(Cot(w),2))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm108 = (C0i(cc0,m12,m22,m32,Mh12,Mh22,MW2))*(1.0*(-(EL*Cos(theta)*Csc(w)*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta))))/(1.0*(64*Power(Pi,2)*vs));
	ComplexType PVTerm109 = (C0i(cc0,m12,m32,m22,MB2,MB2,MT2))*(1.0*(-3*Power(EL,2)*MB2*(m12 + m22 - m32 + 4*MB2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm110 = (C0i(cc0,m12,m32,m22,Mh12,Mh12,MW2))*(1.0*(-3*EL*Csc(w)*(MW*Power(Cos(theta),2)*(d2*Power(vs,2) + 2*(Mh12 - Mh22)*Power(Cos(theta),2))*Power(Sin(theta),3) + EL*Mh12*vs*Csc(w)*(-1 + Power(Cos(theta),7) + Power(Cos(theta),5)*Power(Sin(theta),2)))))/(1.0*(32*Power(Pi,2)*vs));
	ComplexType PVTerm111 = (C0i(cc0,m12,m32,m22,Mh12,Mh22,MW2))*(1.0*(-(EL*Cos(theta)*Csc(w)*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta))))/(1.0*(64*Power(Pi,2)*vs));
	ComplexType PVTerm112 = (C0i(cc0,m12,m32,m22,Mh22,Mh22,MW2))*(1.0*(-(EL*Csc(w)*Power(Sin(theta),2)*(MW*(Mh12 - Mh22 + 6*d2*Power(vs,2))*Cos(theta) + 3*(Mh12 - Mh22)*MW*Cos(3*theta) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta))*Sin(2*theta))))/(1.0*(128*Power(Pi,2)*vs));
	ComplexType PVTerm113 = (C0i(cc0,m12,m32,m22,MW2,MW2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*((2*m12 + m22 - m32 + 2*MW2)*Cot(w)*Csc(w) + 2*(2*m12 - 3*m22 - 2*(m32 + MW2))*Power(Cot(w),3)*Csc(w) + Mh12*Sec(w)))))/(1.0*(32*Power(Pi,2)*Cot(w)*Csc(w)));
	ComplexType PVTerm114 = (C0i(cc0,m22,m12,m32,0,MW2,MW2))*(1.0*(Power(EL,2)*(3*m12 + 2*m22 - 4*m32 + Mh12 - 6*MW2)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm115 = (C0i(cc0,m22,m12,m32,MB2,MT2,MT2))*(1.0*(-3*Power(EL,2)*MT2*(-m22 + MT2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(8*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm116 = (C0i(cc0,m22,m12,m32,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*MW2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm117 = (C0i(cc0,m22,m12,m32,Mh22,MW2,MW2))*(1.0*(-(Power(EL,2)*MW2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm118 = (C0i(cc0,m22,m12,m32,MW2,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),2)*(5*m12 - m22 - 5*m32 - 4*MW2*Power(Sec(w),2) + 2*MW2*Power(Tan(w),4)))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm119 = (C0i(cc00,m12,m22,m32,Mh12,Mh22,MW2))*(1.0*(EL*Cos(theta)*Csc(w)*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta)))/(1.0*(64*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm120 = (C0i(cc00,m12,m32,m22,MB2,MB2,MT2))*(1.0*(-3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(8*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm121 = (C0i(cc00,m12,m32,m22,Mh12,Mh12,MW2))*(1.0*(3*EL*Csc(w)*(MW*Power(Cos(theta),2)*(d2*Power(vs,2) + 2*(Mh12 - Mh22)*Power(Cos(theta),2))*Power(Sin(theta),3) + EL*Mh12*vs*Csc(w)*(-1 + Power(Cos(theta),7) + Power(Cos(theta),5)*Power(Sin(theta),2)))))/(1.0*(32*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm122 = (C0i(cc00,m12,m32,m22,Mh12,Mh22,MW2))*(1.0*(EL*Cos(theta)*Csc(w)*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta)))/(1.0*(64*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm123 = (C0i(cc00,m12,m32,m22,Mh22,Mh22,MW2))*(1.0*(EL*Csc(w)*Power(Sin(theta),2)*(MW*(Mh12 - Mh22 + 6*d2*Power(vs,2))*Cos(theta) + 3*(Mh12 - Mh22)*MW*Cos(3*theta) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta))*Sin(2*theta)))/(1.0*(128*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm124 = (C0i(cc00,m12,m32,m22,MW2,MW2,MZ2))*(1.0*(-(Power(EL,2)*(Mh12 + 10*MW2 + 8*MW2*Cos(2*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm125 = (C0i(cc00,m22,m12,m32,0,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(2*Power(Pi,2)));
	ComplexType PVTerm126 = (C0i(cc00,m22,m12,m32,MB2,MT2,MT2))*(1.0*(-3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(8*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm127 = (C0i(cc00,m22,m12,m32,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12 + 2*MW2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm128 = (C0i(cc00,m22,m12,m32,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12 + 2*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm129 = (C0i(cc00,m22,m12,m32,MW2,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(Mh12 + 20*MW2 + (Mh12 + 16*MW2)*Cos(2*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm130 = (C0i(cc1,m12,m32,m22,MB2,MB2,MT2))*(1.0*(-3*Power(EL,2)*(2*m12 + m22 - m32)*MB2*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm131 = (C0i(cc1,m12,m32,m22,MW2,MW2,MZ2))*(1.0*(-(Power(EL,2)*(m12 + (3*m12 + 2*m22 - 2*m32)*Cos(2*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm132 = (C0i(cc1,m22,m12,m32,0,MW2,MW2))*(1.0*(Power(EL,2)*(m12 + 9*m22 - m32)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm133 = (C0i(cc1,m22,m12,m32,MB2,MT2,MT2))*(1.0*(3*Power(EL,2)*(m12 + 5*m22 - m32)*MT2*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm134 = (C0i(cc1,m22,m12,m32,MW2,MZ2,MZ2))*(1.0*(Power(EL,2)*(m12 + 5*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm135 = (C0i(cc2,m12,m32,m22,MB2,MB2,MT2))*(1.0*(-3*Power(EL,2)*(m12 + 5*m22 - m32)*MB2*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm136 = (C0i(cc2,m12,m32,m22,MW2,MW2,MZ2))*(1.0*(-(Power(EL,2)*(m12 + m22 - m32 + (m12 + 9*m22 - m32)*Cos(2*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm137 = (C0i(cc2,m22,m12,m32,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(5*m12 - 5*m22 - 3*m32)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm138 = (C0i(cc2,m22,m12,m32,MB2,MT2,MT2))*(1.0*(-3*Power(EL,2)*(3*m12 - 3*m22 - m32)*MT2*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm139 = (C0i(cc2,m22,m12,m32,MW2,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(3*m12 - 3*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80 + PVTerm81 + PVTerm82 + PVTerm83 + PVTerm84 + PVTerm85 + PVTerm86 + PVTerm87 + PVTerm88 + PVTerm89 + PVTerm90 + PVTerm91 + PVTerm92 + PVTerm93 + PVTerm94 + PVTerm95 + PVTerm96 + PVTerm97 + PVTerm98 + PVTerm99 + PVTerm100 + PVTerm101 + PVTerm102 + PVTerm103 + PVTerm104 + PVTerm105 + PVTerm106 + PVTerm107 + PVTerm108 + PVTerm109 + PVTerm110 + PVTerm111 + PVTerm112 + PVTerm113 + PVTerm114 + PVTerm115 + PVTerm116 + PVTerm117 + PVTerm118 + PVTerm119 + PVTerm120 + PVTerm121 + PVTerm122 + PVTerm123 + PVTerm124 + PVTerm125 + PVTerm126 + PVTerm127 + PVTerm128 + PVTerm129 + PVTerm130 + PVTerm131 + PVTerm132 + PVTerm133 + PVTerm134 + PVTerm135 + PVTerm136 + PVTerm137 + PVTerm138 + PVTerm139;

}

ComplexType dKappahW_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappahZ(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(-(Power(EL,2)*(-4 - 43*Cos(2*w) + 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(-(Power(EL,2)*(7 + Cos(2*w) + 8*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(Power(EL,2)*(10 + Cos(2*w) + 34*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(7 + 11*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,m12,Mh12,Mh12))*(1.0*(-3*EL*Csc(w)*(MW*Power(Cos(theta),2)*(d2*Power(vs,2) + 2*(Mh12 - Mh22)*Power(Cos(theta),2))*Power(Sin(theta),3) + EL*Mh12*vs*Csc(w)*(-1 + Power(Cos(theta),7) + Power(Cos(theta),5)*Power(Sin(theta),2)))))/(1.0*(128*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm9 = (B0i(bb0,m12,Mh12,Mh22))*(1.0*(EL*Csc(w)*(-2*EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + 2*MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Power(Sin(2*theta),2)))/(1.0*(512*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm10 = (B0i(bb0,m12,Mh22,Mh22))*(1.0*(-(EL*Csc(w)*Power(Sin(theta),2)*(MW*(Mh12 - Mh22 + 6*d2*Power(vs,2))*Cos(theta) + 3*(Mh12 - Mh22)*MW*Cos(3*theta) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta))*Sin(2*theta))))/(1.0*(512*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm11 = (B0i(bb0,m12,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12 + 18*MW2 + 24*MW2*Cos(2*w) + (Mh12 + 6*MW2)*Cos(4*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm12 = (B0i(bb0,m12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm13 = (B0i(bb0,m22,Mh12,MZ2))*(1.0*(Power(EL,2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm14 = (B0i(bb0,m22,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,m22,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*Power(Sin(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,m32,MB2,MB2))*(1.0*(-(Power(EL,2)*MB2*(15 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,m32,Mh12,MZ2))*(1.0*(Power(EL,2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm18 = (B0i(bb0,m32,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm19 = (B0i(bb0,m32,MT2,MT2))*(1.0*(-(Power(EL,2)*MT2*(9 - 2*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(24*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm20 = (B0i(bb0,m32,MW2,MW2))*(1.0*(-(Power(EL,2)*(11 + 20*Cos(2*w) + Cos(4*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm24 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm25 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm26 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Power(Cos(w),2)*(-32*MW2*Power(Csc(2*w),4) + (Cos(theta)*Power(Csc(w),4)*Power(Sec(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(Mh12 - Mh22))))/(1.0*(256*Power(MW,2)*Power(Pi,2)*Power(Csc(w),2)));
	ComplexType PVTerm29 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm31 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm32 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm33 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm34 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*(-1 + 3*Cos(2*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*(-1 + 3*Cos(2*w))*Power(Csc(w),4)*Power(Sec(w),2))))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*(-1 + 3*Cos(2*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta),2))))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(-3*Power(EL,2)*MT2*(-1 + 3*Cos(2*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(27 + 12*Cos(2*w) + 17*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb00,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + 3*Cos(2*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(-6 + 35*Cos(2*w) + 4*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(3 + 2*Cos(theta) + Cos(2*theta))*(-1 + 3*Cos(2*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*(-1 + 3*Cos(2*w))*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(-15 + 37*Cos(2*w) - 10*Cos(4*w) + 6*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(96*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(10 + 35*Cos(2*w) + 18*Cos(4*w) + 9*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm69 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm70 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm71 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm72 = (B0i(bb1,MW2,0,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cot(w),2))*Power(Sin(theta/2.),2)))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Cos(w),2)*(64*Power(EL,2)*Power(Mh12,2)*Power(Csc(2*w),2) - (Cos(theta)*Power(Sec(w),2)*Power(3*EL*Mh12*vs*Cos(theta)*Csc(w) + EL*Mh12*vs*Cos(3*theta)*Csc(w) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm81 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(-2*EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + 9*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Sin(theta) + MW*Sin(theta)*(-Mh12 + Mh22 + 6*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(512*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm82 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(MW*Cos(theta)*(-Mh12 + Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta)) + EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm83 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm84 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm85 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm86 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm87 = (B0i(dbb0,MZ2,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Power(Csc(2*w),2)*Power(Sin(theta/2.),2)))/(1.0*(4*Power(Pi,2)));
	ComplexType PVTerm88 = (B0i(dbb0,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*MW2*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm89 = (B0i(dbb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*MW2*Cos(theta)*Power(Csc(w),2)*Power(Sec(w),4)*Power(Sin(theta),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm90 = (B0i(dbb0,MZ2,MT2,MT2))*(1.0*(-3*Power(EL,2)*MT2*Power(Csc(2*w),2)*Power(Sin(theta/2.),2)))/(1.0*(4*Power(Pi,2)));
	ComplexType PVTerm91 = (B0i(dbb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*MW2*(5 + 9*Cos(2*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm92 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm93 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm94 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm95 = (B0i(dbb00,MZ2,0,0))*(1.0*(3*Power(EL,2)*Power(Csc(2*w),2)*Power(Sin(theta/2.),2)))/(1.0*(2*Power(Pi,2)));
	ComplexType PVTerm96 = (B0i(dbb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm97 = (B0i(dbb00,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm98 = (B0i(dbb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(4*Power(Pi,2)));
	ComplexType PVTerm99 = (B0i(dbb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm100 = (B0i(dbb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm101 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm102 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm103 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm104 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm105 = (B0i(dbb1,MZ2,0,0))*(1.0*(-3*Power(EL,2)*MW2*Power(Csc(2*w),4)*Power(Sin(theta/2.),2)*Power(Sin(w),2)))/(1.0*(Power(Pi,2)));
	ComplexType PVTerm106 = (B0i(dbb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*MW2*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),4)*Power(Sin(theta/2.),2))))/(1.0*(48*Power(Pi,2)));
	ComplexType PVTerm107 = (B0i(dbb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*MW2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),4)*Power(Sin(theta/2.),2))))/(1.0*(48*Power(Pi,2)));
	ComplexType PVTerm108 = (B0i(dbb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*MW2*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm109 = (C0i(cc0,m12,m22,m32,Mh12,Mh22,MZ2))*(1.0*(-(EL*Cos(theta)*Csc(w)*Power(Sec(w),2)*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta))))/(1.0*(64*Power(Pi,2)*vs));
	ComplexType PVTerm110 = (C0i(cc0,m12,m32,m22,Mh12,Mh12,MZ2))*(1.0*(-3*EL*Csc(w)*Power(Sec(w),2)*(MW*Power(Cos(theta),2)*(d2*Power(vs,2) + 2*(Mh12 - Mh22)*Power(Cos(theta),2))*Power(Sin(theta),3) + EL*Mh12*vs*Csc(w)*(-1 + Power(Cos(theta),7) + Power(Cos(theta),5)*Power(Sin(theta),2)))))/(1.0*(32*Power(Pi,2)*vs));
	ComplexType PVTerm111 = (C0i(cc0,m12,m32,m22,Mh12,Mh22,MZ2))*(1.0*(-(EL*Cos(theta)*Csc(w)*Power(Sec(w),2)*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta))))/(1.0*(64*Power(Pi,2)*vs));
	ComplexType PVTerm112 = (C0i(cc0,m12,m32,m22,Mh22,Mh22,MZ2))*(1.0*(-(EL*Csc(w)*Power(Sec(w),2)*Power(Sin(theta),2)*(MW*(Mh12 - Mh22 + 6*d2*Power(vs,2))*Cos(theta) + 3*(Mh12 - Mh22)*MW*Cos(3*theta) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta))*Sin(2*theta))))/(1.0*(128*Power(Pi,2)*vs));
	ComplexType PVTerm113 = (C0i(cc0,m22,m12,m32,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*MW2*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm114 = (C0i(cc0,m22,m12,m32,Mh22,MZ2,MZ2))*(1.0*(-(Power(EL,2)*MW2*Cos(theta)*Power(Csc(w),2)*Power(Sec(w),4)*Power(Sin(theta),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm115 = (C0i(cc0,m22,m32,m12,MB2,MB2,MB2))*(1.0*(Power(EL,2)*MB2*(-1 + Cos(theta))*(6*(m12 + m22 - m32 + 6*MB2) + 2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w))*Power(Csc(w),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm116 = (C0i(cc0,m22,m32,m12,MT2,MT2,MT2))*(1.0*(Power(EL,2)*MT2*(-1 + Cos(theta))*(9*(m12 + m22 - m32 + 4*MT2) - 4*(m12 + m22 - m32)*Cos(2*w) + 4*(m12 + m22 - m32)*Cos(4*w))*Power(Csc(w),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm117 = (C0i(cc0,m22,m32,m12,MW2,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(Mh12 + 2*MW2 - (2*m12 + m22 - m32 + 2*MW2)*Power(Cot(w),2) + (4*m12 - 6*m22 - 4*m32 - 4*MW2)*Power(Cot(w),4))*Power(Sin(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm118 = (C0i(cc00,m12,m22,m32,Mh12,Mh22,MZ2))*(1.0*(EL*Cos(theta)*Csc(w)*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta)))/(1.0*(64*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm119 = (C0i(cc00,m12,m32,m22,Mh12,Mh12,MZ2))*(1.0*(3*EL*Csc(w)*(MW*Power(Cos(theta),2)*(d2*Power(vs,2) + 2*(Mh12 - Mh22)*Power(Cos(theta),2))*Power(Sin(theta),3) + EL*Mh12*vs*Csc(w)*(-1 + Power(Cos(theta),7) + Power(Cos(theta),5)*Power(Sin(theta),2)))))/(1.0*(32*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm120 = (C0i(cc00,m12,m32,m22,Mh12,Mh22,MZ2))*(1.0*(EL*Cos(theta)*Csc(w)*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta)))/(1.0*(64*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm121 = (C0i(cc00,m12,m32,m22,Mh22,Mh22,MZ2))*(1.0*(EL*Csc(w)*Power(Sin(theta),2)*(MW*(Mh12 - Mh22 + 6*d2*Power(vs,2))*Cos(theta) + 3*(Mh12 - Mh22)*MW*Cos(3*theta) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta))*Sin(2*theta)))/(1.0*(128*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm122 = (C0i(cc00,m22,m12,m32,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Mh12 + 2*MW2*Power(Sec(w),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm123 = (C0i(cc00,m22,m12,m32,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12 + 2*MW2*Power(Sec(w),2))*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm124 = (C0i(cc00,m22,m32,m12,MB2,MB2,MB2))*(1.0*(Power(EL,2)*MB2*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm125 = (C0i(cc00,m22,m32,m12,MT2,MT2,MT2))*(1.0*(Power(EL,2)*MT2*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm126 = (C0i(cc00,m22,m32,m12,MW2,MW2,MW2))*(1.0*(-(Power(EL,2)*(Mh12 + 14*MW2 + 16*MW2*Cos(2*w) + (Mh12 + 6*MW2)*Cos(4*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm127 = (C0i(cc1,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-(Power(EL,2)*MB2*(9*m12 + 33*m22 - 9*m32 + 8*m22*Cos(2*w) + 4*m22*Cos(4*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm128 = (C0i(cc1,m22,m32,m12,MT2,MT2,MT2))*(1.0*(-(Power(EL,2)*MT2*(9*(m12 + 5*m22 - m32) - 16*m22*Cos(2*w) + 16*m22*Cos(4*w))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm129 = (C0i(cc1,m22,m32,m12,MW2,MW2,MW2))*(1.0*(-(Power(EL,2)*(m12 + 9*m22 - m32 + (m12 + m22 - m32)*Cos(2*w))*Power(Cot(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm130 = (C0i(cc2,m22,m32,m12,MB2,MB2,MB2))*(1.0*(Power(EL,2)*MB2*(-1 + Cos(theta))*(15*m12 + 6*m22 - 6*m32 + 2*(m12 + m22 - m32)*Cos(2*w) + (m12 + m22 - m32)*Cos(4*w))*Power(Csc(w),2)))/(1.0*(96*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm131 = (C0i(cc2,m22,m32,m12,MT2,MT2,MT2))*(1.0*(Power(EL,2)*MT2*(-1 + Cos(theta))*(9*(2*m12 + m22 - m32) - 4*(m12 + m22 - m32)*Cos(2*w) + 4*(m12 + m22 - m32)*Cos(4*w))*Power(Csc(w),2)))/(1.0*(96*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm132 = (C0i(cc2,m22,m32,m12,MW2,MW2,MW2))*(1.0*(-(Power(EL,2)*(3*m12 + 2*m22 - 2*m32 + m12*Cos(2*w))*Power(Cot(w),2)*Power(Sin(theta/2.),2))))/(1.0*(8*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80 + PVTerm81 + PVTerm82 + PVTerm83 + PVTerm84 + PVTerm85 + PVTerm86 + PVTerm87 + PVTerm88 + PVTerm89 + PVTerm90 + PVTerm91 + PVTerm92 + PVTerm93 + PVTerm94 + PVTerm95 + PVTerm96 + PVTerm97 + PVTerm98 + PVTerm99 + PVTerm100 + PVTerm101 + PVTerm102 + PVTerm103 + PVTerm104 + PVTerm105 + PVTerm106 + PVTerm107 + PVTerm108 + PVTerm109 + PVTerm110 + PVTerm111 + PVTerm112 + PVTerm113 + PVTerm114 + PVTerm115 + PVTerm116 + PVTerm117 + PVTerm118 + PVTerm119 + PVTerm120 + PVTerm121 + PVTerm122 + PVTerm123 + PVTerm124 + PVTerm125 + PVTerm126 + PVTerm127 + PVTerm128 + PVTerm129 + PVTerm130 + PVTerm131 + PVTerm132;

}

ComplexType dKappahZ_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappahb(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(-(Power(EL,2)*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(5 + 13*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,m32,0,MB2))*(1.0*(-((Power(EL,2) + 48*Alfas*Pi)*Power(Sin(theta/2.),2))))/(1.0*(18*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,m32,MB2,Mh12))*(1.0*(-(Power(EL,2)*MB2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm10 = (B0i(bb0,m32,MB2,Mh22))*(1.0*(-(Power(EL,2)*MB2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm11 = (B0i(bb0,m32,MB2,MZ2))*(1.0*(-(Power(EL,2)*(9*MB2 + 6*MW2 + (9*MB2 + 8*MW2)*Cos(2*w) + 4*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(576*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm12 = (B0i(bb0,m32,MT2,MW2))*(1.0*(-(Power(EL,2)*(MT2 + MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm13 = (B0i(bb0,MB2,0,MB2))*(1.0*((Power(EL,2) + 48*Alfas*Pi)*Power(Sin(theta/2.),2)))/(1.0*(18*Power(Pi,2)));
	ComplexType PVTerm14 = (B0i(bb0,MB2,MB2,Mh12))*(1.0*(Power(EL,2)*MB2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,MB2,MB2,Mh22))*(1.0*(Power(EL,2)*MB2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,MB2,MB2,MZ2))*(1.0*(Power(EL,2)*(9*MB2 - 12*MW2 + (9*MB2 + 8*MW2)*Cos(2*w) + 4*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(576*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,MB2,MT2,MW2))*(1.0*(-(Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Csc(w),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm18 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm19 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm20 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm24 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm25 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(2*MW2*(-4*Mh12 + 4*Mh22 + Mh12*Cos(theta) - Mh12*Cos(3*theta))*Power(Csc(2*w),2) + Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm26 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm28 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm29 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm30 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm31 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm32 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm33 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm34 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm69 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm70 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm71 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm72 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(dbb0,MB2,0,MB2))*(1.0*(-(MB2*(Power(EL,2) + 48*Alfas*Pi)*Power(Sin(theta/2.),2))))/(1.0*(18*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(dbb0,MB2,MB2,Mh12))*(1.0*(-(Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(dbb0,MB2,MB2,Mh22))*(1.0*(-(Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(dbb0,MB2,MB2,MZ2))*(1.0*(-(Power(EL,2)*MB2*(9*MB2 - 12*MW2 + (9*MB2 + 8*MW2)*Cos(2*w) + 4*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(288*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb0,MB2,MT2,MW2))*(1.0*(Power(EL,2)*MB2*MT2*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm81 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm82 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm83 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm84 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm85 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm86 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm87 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm88 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm89 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm90 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm91 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm92 = (B0i(dbb1,MB2,0,MB2))*(1.0*(MB2*(Power(EL,2) + 48*Alfas*Pi)*Power(Sin(theta/2.),2)))/(1.0*(18*Power(Pi,2)));
	ComplexType PVTerm93 = (B0i(dbb1,MB2,MB2,Mh12))*(1.0*(Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm94 = (B0i(dbb1,MB2,MB2,Mh22))*(1.0*(Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm95 = (B0i(dbb1,MB2,MB2,MZ2))*(1.0*(-(Power(EL,2)*MB2*(9*MB2 + 12*MW2 + (9*MB2 + 4*MW2)*Cos(2*w) + 2*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(288*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm96 = (B0i(dbb1,MB2,MT2,MW2))*(1.0*(-(Power(EL,2)*MB2*(MB2 + MT2 + 2*MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm97 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm98 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm99 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm100 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm101 = (C0i(cc0,m12,m32,m22,MB2,MB2,Mh12))*(1.0*(-(Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm102 = (C0i(cc0,m12,m32,m22,MB2,MB2,Mh22))*(1.0*(-(Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm103 = (C0i(cc0,m12,m32,m22,MB2,MB2,MZ2))*(1.0*(-(Power(EL,2)*MB2*(9*MB2 - 12*MW2 + (9*MB2 + 8*MW2)*Cos(2*w) + 4*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(288*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm104 = (C0i(cc0,m12,m32,m22,MT2,MT2,MW2))*(1.0*(Power(EL,2)*Power(MT2,2)*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm105 = (C0i(cc0,m22,m12,m32,0,MB2,MB2))*(1.0*(-((m12 + m22 - m32 - 4*MB2)*(Power(EL,2) + 48*Alfas*Pi)*(-1 + Cos(theta)))))/(1.0*(72*Power(Pi,2)));
	ComplexType PVTerm106 = (C0i(cc0,m22,m12,m32,MB2,Mh12,Mh12))*(1.0*(-3*EL*MB2*Csc(w)*(MW*Power(Cos(theta),2)*(d2*Power(vs,2) + 2*(Mh12 - Mh22)*Power(Cos(theta),2))*Power(Sin(theta),3) + EL*Mh12*vs*Csc(w)*(-1 + Power(Cos(theta),7) + Power(Cos(theta),5)*Power(Sin(theta),2)))))/(1.0*(64*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm107 = (C0i(cc0,m22,m12,m32,MB2,Mh12,Mh22))*(1.0*(EL*MB2*Cos(theta)*Csc(w)*Sin(theta)*(-2*EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + 2*MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta)))/(1.0*(256*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm108 = (C0i(cc0,m22,m12,m32,MB2,Mh22,Mh22))*(1.0*(-(EL*MB2*Csc(w)*Power(Sin(theta),2)*(MW*(Mh12 - Mh22 + 6*d2*Power(vs,2))*Cos(theta) + 3*(Mh12 - Mh22)*MW*Cos(3*theta) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta))*Sin(2*theta))))/(1.0*(256*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm109 = (C0i(cc0,m22,m12,m32,MB2,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-3*(-9*MB2*Mh12 + 8*Power(MW,4) + 12*m22*MW2) + 4*(9*MB2*Mh12 + 16*Power(MW,4) - 9*m22*MW2)*Cos(2*w) + (9*MB2*Mh12 + 32*Power(MW,4))*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),4)*Power(Sin(theta/2.),2))))/(1.0*(2304*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm110 = (C0i(cc0,m22,m12,m32,MT2,MW2,MW2))*(1.0*(-(Power(EL,2)*(Mh12*MT2 + Power(MW,4) - m22*MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm111 = (C0i(cc0,m32,m12,m22,MB2,Mh12,Mh22))*(1.0*(EL*MB2*Cos(theta)*Csc(w)*Sin(theta)*(-2*EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + 2*MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta)))/(1.0*(256*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm112 = (C0i(cc1,m12,m32,m22,MB2,MB2,Mh12))*(1.0*(-(Power(EL,2)*m12*MB2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm113 = (C0i(cc1,m12,m32,m22,MB2,MB2,Mh22))*(1.0*(-(Power(EL,2)*m12*MB2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm114 = (C0i(cc1,m12,m32,m22,MB2,MB2,MZ2))*(1.0*(-(Power(EL,2)*m12*(9*MB2 - 12*MW2 + (9*MB2 + 8*MW2)*Cos(2*w) + 4*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(576*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm115 = (C0i(cc1,m12,m32,m22,MT2,MT2,MW2))*(1.0*(Power(EL,2)*m12*MT2*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm116 = (C0i(cc1,m22,m12,m32,0,MB2,MB2))*(1.0*((m12 + m22 - m32)*(Power(EL,2) + 48*Alfas*Pi)*Power(Sin(theta/2.),2)))/(1.0*(36*Power(Pi,2)));
	ComplexType PVTerm117 = (C0i(cc1,m22,m12,m32,MB2,MZ2,MZ2))*(1.0*(Power(EL,2)*(m12 + 5*m22 - m32)*Power(Csc(2*w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm118 = (C0i(cc1,m22,m12,m32,MT2,MW2,MW2))*(1.0*(Power(EL,2)*(m12 + 5*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm119 = (C0i(cc2,m12,m32,m22,MB2,MB2,Mh12))*(1.0*(Power(EL,2)*(m12 + m22 - m32)*MB2*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm120 = (C0i(cc2,m12,m32,m22,MB2,MB2,Mh22))*(1.0*(-(Power(EL,2)*(m12 + m22 - m32)*MB2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm121 = (C0i(cc2,m12,m32,m22,MB2,MB2,MZ2))*(1.0*(Power(EL,2)*(m12 + m22 - m32)*(-1 + Cos(theta))*(9*MB2 - 12*MW2 + (9*MB2 + 8*MW2)*Cos(2*w) + 4*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)))/(1.0*(2304*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm122 = (C0i(cc2,m12,m32,m22,MT2,MT2,MW2))*(1.0*(-(Power(EL,2)*(m12 + m22 - m32)*MT2*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm123 = (C0i(cc2,m22,m12,m32,0,MB2,MB2))*(1.0*(-((m12 - m22 + m32)*(Power(EL,2) + 48*Alfas*Pi)*Power(Sin(theta/2.),2))))/(1.0*(36*Power(Pi,2)));
	ComplexType PVTerm124 = (C0i(cc2,m22,m12,m32,MB2,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(3*m12 - 3*m22 - m32)*Power(Csc(2*w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm125 = (C0i(cc2,m22,m12,m32,MT2,MW2,MW2))*(1.0*(-(Power(EL,2)*(3*m12 - 3*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80 + PVTerm81 + PVTerm82 + PVTerm83 + PVTerm84 + PVTerm85 + PVTerm86 + PVTerm87 + PVTerm88 + PVTerm89 + PVTerm90 + PVTerm91 + PVTerm92 + PVTerm93 + PVTerm94 + PVTerm95 + PVTerm96 + PVTerm97 + PVTerm98 + PVTerm99 + PVTerm100 + PVTerm101 + PVTerm102 + PVTerm103 + PVTerm104 + PVTerm105 + PVTerm106 + PVTerm107 + PVTerm108 + PVTerm109 + PVTerm110 + PVTerm111 + PVTerm112 + PVTerm113 + PVTerm114 + PVTerm115 + PVTerm116 + PVTerm117 + PVTerm118 + PVTerm119 + PVTerm120 + PVTerm121 + PVTerm122 + PVTerm123 + PVTerm124 + PVTerm125;

}

ComplexType dKappahb_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappahc(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(-(Power(EL,2)*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(5 + 13*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm10 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm11 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm12 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm13 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm14 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(2*MW2*(-4*Mh12 + 4*Mh22 + Mh12*Cos(theta) - Mh12*Cos(3*theta))*Power(Csc(2*w),2) + Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm18 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm19 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm20 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm24 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm25 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm26 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm29 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm31 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm32 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm33 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm34 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm69 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm70 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm71 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm72 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80;

}

ComplexType dKappahc_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappahd(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(-(Power(EL,2)*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(5 + 13*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm10 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm11 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm12 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm13 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm14 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(2*MW2*(-4*Mh12 + 4*Mh22 + Mh12*Cos(theta) - Mh12*Cos(3*theta))*Power(Csc(2*w),2) + Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm18 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm19 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm20 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm24 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm25 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm26 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm29 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm31 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm32 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm33 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm34 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm69 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm70 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm71 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm72 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80;

}

ComplexType dKappahd_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappahele(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(-(Power(EL,2)*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(5 + 13*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,m32,0,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm10 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm11 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm12 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm13 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm14 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(2*MW2*(-4*Mh12 + 4*Mh22 + Mh12*Cos(theta) - Mh12*Cos(3*theta))*Power(Csc(2*w),2) + Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm18 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm19 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm20 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm24 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm25 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm26 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm29 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm31 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm32 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm33 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm34 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm69 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm70 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm71 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm72 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm81 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm82 = (C0i(cc0,m22,m12,m32,0,MW2,MW2))*(1.0*(Power(EL,2)*(m22 - MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm83 = (C0i(cc1,m22,m12,m32,0,MW2,MW2))*(1.0*(Power(EL,2)*(m12 + 5*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm84 = (C0i(cc2,m22,m12,m32,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(3*m12 - 3*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80 + PVTerm81 + PVTerm82 + PVTerm83 + PVTerm84;

}

ComplexType dKappahele_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappahmuon(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(-(Power(EL,2)*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(5 + 13*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,m32,0,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm10 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm11 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm12 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm13 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm14 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(2*MW2*(-4*Mh12 + 4*Mh22 + Mh12*Cos(theta) - Mh12*Cos(3*theta))*Power(Csc(2*w),2) + Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm18 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm19 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm20 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm24 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm25 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm26 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm29 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm31 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm32 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm33 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm34 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm69 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm70 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm71 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm72 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm81 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm82 = (C0i(cc0,m22,m12,m32,0,MW2,MW2))*(1.0*(Power(EL,2)*(m22 - MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm83 = (C0i(cc1,m22,m12,m32,0,MW2,MW2))*(1.0*(Power(EL,2)*(m12 + 5*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm84 = (C0i(cc2,m22,m12,m32,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(3*m12 - 3*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80 + PVTerm81 + PVTerm82 + PVTerm83 + PVTerm84;

}

ComplexType dKappahmuon_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappahs(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(-(Power(EL,2)*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(5 + 13*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm10 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm11 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm12 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm13 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm14 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(2*MW2*(-4*Mh12 + 4*Mh22 + Mh12*Cos(theta) - Mh12*Cos(3*theta))*Power(Csc(2*w),2) + Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm18 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm19 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm20 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm24 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm25 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm26 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm29 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm31 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm32 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm33 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm34 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm69 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm70 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm71 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm72 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80;

}

ComplexType dKappahs_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappaht(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(-(Power(EL,2)*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(5 + 13*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,m32,0,MT2))*(1.0*(-2*(Power(EL,2) + 12*Alfas*Pi)*Power(Sin(theta/2.),2)))/(1.0*(9*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,m32,MB2,MW2))*(1.0*(-(Power(EL,2)*(MB2 + MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm10 = (B0i(bb0,m32,Mh12,MT2))*(1.0*(-(Power(EL,2)*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm11 = (B0i(bb0,m32,Mh22,MT2))*(1.0*(-(Power(EL,2)*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm12 = (B0i(bb0,m32,MT2,MZ2))*(1.0*(-(Power(EL,2)*(9*(MT2 + 2*MW2) + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(576*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm13 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm14 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm15 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm16 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm17 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm18 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm19 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm20 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(2*MW2*(-4*Mh12 + 4*Mh22 + Mh12*Cos(theta) - Mh12*Cos(3*theta))*Power(Csc(2*w),2) + Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm24 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm25 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm26 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm29 = (B0i(bb0,MT2,0,MT2))*(1.0*(2*(Power(EL,2) + 12*Alfas*Pi)*Power(Sin(theta/2.),2)))/(1.0*(9*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,MT2,MB2,MW2))*(1.0*(-(Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Csc(w),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm31 = (B0i(bb0,MT2,Mh12,MT2))*(1.0*(Power(EL,2)*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm32 = (B0i(bb0,MT2,Mh22,MT2))*(1.0*(Power(EL,2)*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm33 = (B0i(bb0,MT2,MT2,MZ2))*(1.0*(Power(EL,2)*(9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(576*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm34 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm69 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm70 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm71 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm72 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm79 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm80 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm81 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm82 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm83 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm84 = (B0i(dbb0,MT2,0,MT2))*(1.0*(-2*MT2*(Power(EL,2) + 12*Alfas*Pi)*Power(Sin(theta/2.),2)))/(1.0*(9*Power(Pi,2)));
	ComplexType PVTerm85 = (B0i(dbb0,MT2,MB2,MW2))*(1.0*(Power(EL,2)*MB2*MT2*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm86 = (B0i(dbb0,MT2,Mh12,MT2))*(1.0*(-(Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm87 = (B0i(dbb0,MT2,Mh22,MT2))*(1.0*(-(Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm88 = (B0i(dbb0,MT2,MT2,MZ2))*(1.0*(-(Power(EL,2)*MT2*(9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(288*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm89 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm90 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm91 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm92 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm93 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm94 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm95 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm96 = (B0i(dbb1,MT2,0,MT2))*(1.0*(2*MT2*(Power(EL,2) + 12*Alfas*Pi)*Power(Sin(theta/2.),2)))/(1.0*(9*Power(Pi,2)));
	ComplexType PVTerm97 = (B0i(dbb1,MT2,MB2,MW2))*(1.0*(-(Power(EL,2)*MT2*(MB2 + MT2 + 2*MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm98 = (B0i(dbb1,MT2,Mh12,MT2))*(1.0*(-(Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm99 = (B0i(dbb1,MT2,Mh22,MT2))*(1.0*(-(Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm100 = (B0i(dbb1,MT2,MT2,MZ2))*(1.0*(-(Power(EL,2)*MT2*(9*(MT2 + 2*MW2) + (9*MT2 - 8*MW2)*Cos(2*w) + 8*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(288*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm101 = (C0i(cc0,m12,m22,m32,Mh12,Mh22,MT2))*(1.0*(EL*MT2*Cos(theta)*Csc(w)*Sin(theta)*(-2*EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + 2*MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta)))/(1.0*(256*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm102 = (C0i(cc0,m12,m32,m22,MB2,MB2,MW2))*(1.0*(Power(EL,2)*Power(MB2,2)*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm103 = (C0i(cc0,m12,m32,m22,Mh12,Mh12,MT2))*(1.0*(-3*EL*MT2*Csc(w)*(MW*Power(Cos(theta),2)*(d2*Power(vs,2) + 2*(Mh12 - Mh22)*Power(Cos(theta),2))*Power(Sin(theta),3) + EL*Mh12*vs*Csc(w)*(-1 + Power(Cos(theta),7) + Power(Cos(theta),5)*Power(Sin(theta),2)))))/(1.0*(64*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm104 = (C0i(cc0,m12,m32,m22,Mh12,Mh22,MT2))*(1.0*(EL*MT2*Cos(theta)*Csc(w)*Sin(theta)*(-2*EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + 2*MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*Sin(2*theta)))/(1.0*(256*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm105 = (C0i(cc0,m12,m32,m22,Mh22,Mh22,MT2))*(1.0*(-(EL*MT2*Csc(w)*Power(Sin(theta),2)*(MW*(Mh12 - Mh22 + 6*d2*Power(vs,2))*Cos(theta) + 3*(Mh12 - Mh22)*MW*Cos(3*theta) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta))*Sin(2*theta))))/(1.0*(256*Power(MW,2)*Power(Pi,2)*vs));
	ComplexType PVTerm106 = (C0i(cc0,m12,m32,m22,MT2,MT2,MZ2))*(1.0*(-(Power(EL,2)*MT2*(9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(288*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm107 = (C0i(cc0,m22,m12,m32,0,MT2,MT2))*(1.0*(-((m12 + m22 - m32 - 4*MT2)*(Power(EL,2) + 12*Alfas*Pi)*(-1 + Cos(theta)))))/(1.0*(18*Power(Pi,2)));
	ComplexType PVTerm108 = (C0i(cc0,m22,m12,m32,MB2,MW2,MW2))*(1.0*(-(Power(EL,2)*(MB2*Mh12 + Power(MW,4) - m22*MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm109 = (C0i(cc0,m22,m12,m32,Mh12,MT2,MT2))*(1.0*(Power(EL,2)*MT2*(-m12 - m22 + m32 + 4*MT2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm110 = (C0i(cc0,m22,m12,m32,Mh22,MT2,MT2))*(1.0*(Power(EL,2)*(m12 + m22 - m32 - 4*MT2)*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm111 = (C0i(cc0,m22,m12,m32,MT2,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(9*(3*Mh12*MT2 + 8*Power(MW,4) - 4*m22*MW2) + 4*(9*Mh12*MT2 - 32*Power(MW,4) - 9*m22*MW2)*Cos(2*w) + (9*Mh12*MT2 + 128*Power(MW,4))*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),4)*Power(Sin(theta/2.),2))))/(1.0*(2304*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm112 = (C0i(cc1,m12,m32,m22,MB2,MB2,MW2))*(1.0*(Power(EL,2)*m12*MB2*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm113 = (C0i(cc1,m12,m32,m22,MT2,MT2,MZ2))*(1.0*(-(Power(EL,2)*m12*(9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(576*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm114 = (C0i(cc1,m22,m12,m32,0,MT2,MT2))*(1.0*(-((m12 + m22 - m32)*(Power(EL,2) + 12*Alfas*Pi)*(-1 + Cos(theta)))))/(1.0*(18*Power(Pi,2)));
	ComplexType PVTerm115 = (C0i(cc1,m22,m12,m32,MB2,MW2,MW2))*(1.0*(Power(EL,2)*(m12 + 5*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm116 = (C0i(cc1,m22,m12,m32,Mh12,MT2,MT2))*(1.0*(-(Power(EL,2)*(m12 + m22 - m32)*MT2*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm117 = (C0i(cc1,m22,m12,m32,Mh22,MT2,MT2))*(1.0*(Power(EL,2)*(m12 + m22 - m32)*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm118 = (C0i(cc1,m22,m12,m32,MT2,MZ2,MZ2))*(1.0*(Power(EL,2)*(m12 + 5*m22 - m32)*Power(Csc(2*w),2)*Power(Sin(theta/2.),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm119 = (C0i(cc2,m12,m32,m22,MB2,MB2,MW2))*(1.0*(-(Power(EL,2)*(m12 + m22 - m32)*MB2*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm120 = (C0i(cc2,m12,m32,m22,MT2,MT2,MZ2))*(1.0*(Power(EL,2)*(m12 + m22 - m32)*(-1 + Cos(theta))*(9*MT2 + (9*MT2 - 16*MW2)*Cos(2*w) + 16*MW2*Cos(4*w))*Power(Csc(w),2)*Power(Sec(w),2)))/(1.0*(2304*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm121 = (C0i(cc2,m22,m12,m32,0,MT2,MT2))*(1.0*((m12 - m22 + m32)*(Power(EL,2) + 12*Alfas*Pi)*(-1 + Cos(theta))))/(1.0*(18*Power(Pi,2)));
	ComplexType PVTerm122 = (C0i(cc2,m22,m12,m32,MB2,MW2,MW2))*(1.0*(-(Power(EL,2)*(3*m12 - 3*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm123 = (C0i(cc2,m22,m12,m32,Mh12,MT2,MT2))*(1.0*(Power(EL,2)*(m12 - m22 + m32)*MT2*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm124 = (C0i(cc2,m22,m12,m32,Mh22,MT2,MT2))*(1.0*(-(Power(EL,2)*(m12 - m22 + m32)*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm125 = (C0i(cc2,m22,m12,m32,MT2,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(3*m12 - 3*m22 - m32)*Power(Csc(2*w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80 + PVTerm81 + PVTerm82 + PVTerm83 + PVTerm84 + PVTerm85 + PVTerm86 + PVTerm87 + PVTerm88 + PVTerm89 + PVTerm90 + PVTerm91 + PVTerm92 + PVTerm93 + PVTerm94 + PVTerm95 + PVTerm96 + PVTerm97 + PVTerm98 + PVTerm99 + PVTerm100 + PVTerm101 + PVTerm102 + PVTerm103 + PVTerm104 + PVTerm105 + PVTerm106 + PVTerm107 + PVTerm108 + PVTerm109 + PVTerm110 + PVTerm111 + PVTerm112 + PVTerm113 + PVTerm114 + PVTerm115 + PVTerm116 + PVTerm117 + PVTerm118 + PVTerm119 + PVTerm120 + PVTerm121 + PVTerm122 + PVTerm123 + PVTerm124 + PVTerm125;

}

ComplexType dKappaht_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappahtau(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(-(Power(EL,2)*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(5 + 13*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,m32,0,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm10 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm11 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm12 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm13 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm14 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(2*MW2*(-4*Mh12 + 4*Mh22 + Mh12*Cos(theta) - Mh12*Cos(3*theta))*Power(Csc(2*w),2) + Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm18 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm19 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm20 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm24 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm25 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm26 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm29 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm31 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm32 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm33 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm34 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm69 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm70 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm71 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm72 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm81 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm82 = (C0i(cc0,m22,m12,m32,0,MW2,MW2))*(1.0*(Power(EL,2)*(m22 - MW2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm83 = (C0i(cc1,m22,m12,m32,0,MW2,MW2))*(1.0*(Power(EL,2)*(m12 + 5*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm84 = (C0i(cc2,m22,m12,m32,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(3*m12 - 3*m22 - m32)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80 + PVTerm81 + PVTerm82 + PVTerm83 + PVTerm84;

}

ComplexType dKappahtau_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dKappahu(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = -2*Power(Sin(theta/2.),2);
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(Power(EL,2)*(-16 + 61*Cos(2*w) - 10*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(384*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,Mh12))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,Mh22))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm4 = (A0i(aa0,MT2))*(1.0*(Power(EL,2)*(-29 + 29*Cos(2*w) - 20*Cos(4*w) + 2*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(192*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm5 = (A0i(aa0,MW2))*(1.0*(-(Power(EL,2)*(-22 + Cos(2*w) - 30*Cos(4*w) + 3*Cos(6*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm6 = (A0i(aa0,MZ2))*(1.0*(-(Power(EL,2)*(5 + 13*Cos(2*w) + 6*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,0,MW2,MW2))*(1.0*(-9*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,Mh12,Mh12,Mh12))*(1.0*(3*Sin(theta)*(-(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w)) + MW*(Mh12 - Mh22 + 3*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm10 = (B0i(bb0,Mh12,Mh12,Mh22))*(1.0*(-(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta)))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm11 = (B0i(bb0,Mh12,Mh22,Mh22))*(1.0*(3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm12 = (B0i(bb0,Mh12,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(d2*Power(vs,2) - 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm13 = (B0i(bb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm14 = (B0i(bb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),3) - Mh12*Mh22*Cos(theta)*Power(Sin(theta),2) + MW2*(-2*Mh12 + 2*Mh22 + (Mh12 - 6*MW2)*Sin(theta)*Sin(2*theta)))))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm15 = (B0i(bb0,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*(2*MW2*(-4*Mh12 + 4*Mh22 + Mh12*Cos(theta) - Mh12*Cos(3*theta))*Power(Csc(2*w),2) + Cos(theta)*Power(Csc(w),2)*(2*(Mh12 - Mh22)*MW2*Power(Cos(theta),2)*Power(Sec(w),2) - (Mh12*Mh22 + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm16 = (B0i(bb0,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*Power(MB2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm17 = (B0i(bb0,Mh22,Mh12,Mh12))*(1.0*(3*Sin(theta)*(EL*(2*Mh12 + Mh22)*vs*Cos(theta)*Csc(w) + MW*(-Mh12 + Mh22 - 3*d2*Power(vs,2) - 3*(Mh12 - Mh22)*Cos(2*theta))*Sin(theta))*(EL*Mh12*vs*Power(Cos(theta),5)*Csc(w) + EL*Mh12*vs*Power(Cos(theta),3)*Csc(w)*Power(Sin(theta),2) + d2*MW*Power(vs,2)*Power(Sin(theta),3) + 2*(Mh12 - Mh22)*MW*Power(Cos(theta),2)*Power(Sin(theta),3))*Sin(2*theta)))/(1.0*(512*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm18 = (B0i(bb0,Mh22,Mh12,Mh22))*(1.0*(Sin(theta)*Power(Sin(2*theta),2)*(EL*MW*vs*(2*(2*Power(Mh12,2) - 4*Mh12*Mh22 + 2*Power(Mh22,2) + 9*d2*Mh12*Power(vs,2) + 9*d2*Mh22*Power(vs,2))*Cos(2*theta) + 3*(Mh12 - Mh22)*(Mh12 + Mh22 + 2*d2*Power(vs,2) + 3*(Mh12 + Mh22)*Cos(4*theta)))*Csc(w) - MW2*(7*Power(Mh12,2) - 14*Mh12*Mh22 + 7*Power(Mh22,2) + 18*Power(d2,2)*Power(vs,4) + 36*d2*(Mh12 - Mh22)*Power(vs,2)*Cos(2*theta) + 9*Power(Mh12 - Mh22,2)*Cos(4*theta))*Sin(2*theta) + 2*Power(EL,2)*(2*Power(Mh12,2) + 5*Mh12*Mh22 + 2*Power(Mh22,2))*Power(vs,2)*Power(Csc(w),2)*Sin(2*theta))))/(1.0*(2048*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm19 = (B0i(bb0,Mh22,Mh22,Mh22))*(1.0*(-3*Power(Csc(w),2)*Sin(theta)*Sin(2*theta)*(-12*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + 2*EL*MW*vs*(-((Mh12 - Mh22)*(3*Mh12 - Mh22 - 4*d2*Power(vs,2))*Sin(2*theta)) + 2*(-4*Mh12*Mh22 + 4*Power(Mh22,2) + d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2))*Sin(4*theta) + (Power(Mh12,2) + 4*Mh12*Mh22 - 5*Power(Mh22,2))*Sin(6*theta))*Sin(w) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) + 36*Power(d2,2)*Power(vs,4) + 4*(Mh12 - Mh22)*(Mh12 - Mh22 + 3*d2*Power(vs,2))*Cos(6*theta) + 3*Power(Mh12 - Mh22,2)*Cos(8*theta))*Power(Sin(w),2) - 4*(Cos(2*theta)*(-4*Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) - 2*Mh12*Mh22 + Power(Mh22,2) - 5*d2*Mh12*Power(vs,2) + 5*d2*Mh22*Power(vs,2) - 12*Power(d2,2)*Power(vs,4))*Power(Sin(w),2)) + Cos(4*theta)*(Power(EL,2)*Mh22*(Mh12 + 2*Mh22)*Power(vs,2) + MW2*(Power(Mh12,2) + Power(Mh22,2) + 8*d2*Mh22*Power(vs,2) - 3*Power(d2,2)*Power(vs,4) - 2*Mh12*(Mh22 + 4*d2*Power(vs,2)))*Power(Sin(w),2))))))/(1.0*(16384*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm20 = (B0i(bb0,Mh22,MHA2,MHA2))*(1.0*(Cos(theta)*(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta))*Power(Sin(theta),2)*(-(d2*Power(vs,2)) + 2*(Mh12 - Mh22)*Power(Sin(theta),2))))/(1.0*(256*(Mh12 - Mh22)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm21 = (B0i(bb0,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Power(MT2,2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm22 = (B0i(bb0,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*(Mh12*Mh22 + 12*Power(MW,4) - 2*Mh22*MW2)*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(128*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm23 = (B0i(bb0,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Csc(w),2)*(Mh12*Mh22 - 2*Mh22*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4))*Power(Sin(theta),2)))/(1.0*(256*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm24 = (B0i(bb0,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm25 = (B0i(bb0,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(MT2 - MW2)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm26 = (B0i(bb0,MW2,Mh12,MW2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm27 = (B0i(bb0,MW2,Mh22,MW2))*(1.0*(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm28 = (B0i(bb0,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(12 + 15*Cos(2*w) + 12*Cos(4*w) + Cos(6*w))*Power(Csc(w),4)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(Pi,2)));
	ComplexType PVTerm29 = (B0i(bb0,MZ2,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm30 = (B0i(bb0,MZ2,Mh12,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),4))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm31 = (B0i(bb0,MZ2,Mh22,MZ2))*(1.0*(-(Power(EL,2)*Cos(theta)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm32 = (B0i(bb0,MZ2,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(64*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm33 = (B0i(bb0,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(5 + 9*Cos(2*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm34 = (B0i(bb00,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(1 + 2*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(12*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm35 = (B0i(bb00,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + 4*Cos(2*w))*Power(Sin(theta/2.),2))))/(1.0*(6*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm36 = (B0i(bb00,0,MW2,MW2))*(1.0*(Power(EL,2)*(2 + 3*Cos(2*w))*Power(Sin(theta/2.),2)))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm37 = (B0i(bb00,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(4*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm38 = (B0i(bb00,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm39 = (B0i(bb00,MW2,Mh12,MW2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Cos(2*w)*Power(Csc(w),4))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm40 = (B0i(bb00,MW2,Mh22,MW2))*(1.0*(-(Power(EL,2)*Cos(theta)*Cos(2*w)*Power(Csc(w),4)*Power(Sin(theta),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm41 = (B0i(bb00,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(2 + 5*Cos(2*w) + 2*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm42 = (B0i(bb00,MZ2,0,0))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm43 = (B0i(bb00,MZ2,MB2,MB2))*(1.0*(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm44 = (B0i(bb00,MZ2,Mh12,MZ2))*(1.0*(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Cot(w),2)*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm45 = (B0i(bb00,MZ2,Mh22,MZ2))*(1.0*(Power(EL,2)*Cos(theta)*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm46 = (B0i(bb00,MZ2,MT2,MT2))*(1.0*(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(48*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm47 = (B0i(bb00,MZ2,MW2,MW2))*(1.0*(-(Power(EL,2)*(7 + 8*Cos(2*w) + 3*Cos(4*w))*Power(Cot(w),2)*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm48 = (B0i(bb1,0,MB2,MB2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(24*Power(Pi,2)));
	ComplexType PVTerm49 = (B0i(bb1,0,MT2,MT2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(6*Power(Pi,2)));
	ComplexType PVTerm50 = (B0i(bb1,0,MW2,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm51 = (B0i(bb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm52 = (B0i(bb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*MT2*(-2*Mh12 + 2*Mh22 - Cos(theta)*(-2*Mh12 + Mh22 + Mh22*Cos(2*theta)))*Power(Csc(w),2)))/(1.0*(64*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm53 = (B0i(bb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2))))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm54 = (B0i(bb1,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-2*Mh12 + 3*Mh22 + 2*Mh22*Cos(theta) + Mh22*Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2))))/(1.0*(64*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm55 = (B0i(bb1,Mh22,MB2,MB2))*(1.0*(-3*Power(EL,2)*MB2*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm56 = (B0i(bb1,Mh22,MT2,MT2))*(1.0*(-3*Power(EL,2)*Mh22*MT2*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm57 = (B0i(bb1,Mh22,MW2,MW2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(w),2)*Power(Sin(theta),2)))/(1.0*(32*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm58 = (B0i(bb1,Mh22,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh22*Cos(theta)*Power(Csc(2*w),2)*Power(Sin(theta),2)))/(1.0*(16*(Mh12 - Mh22)*Power(Pi,2)));
	ComplexType PVTerm59 = (B0i(bb1,MW2,0,MW2))*(1.0*(-(Power(EL,2)*(-1 + Cos(theta))*(-1 + Power(Cot(w),2)))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm60 = (B0i(bb1,MW2,MB2,MT2))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Cos(2*w)*Power(Csc(w),4)))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm61 = (B0i(bb1,MW2,MW2,MZ2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),2)*(-1 + Power(Cot(w),2))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm62 = (B0i(bb1,MZ2,0,0))*(1.0*(3*Power(EL,2)*(-1 + Cos(theta))*Power(Csc(w),4)))/(1.0*(64*Power(Pi,2)));
	ComplexType PVTerm63 = (B0i(bb1,MZ2,MB2,MB2))*(1.0*(-(Power(EL,2)*(6 + 2*Cos(2*w) + Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm64 = (B0i(bb1,MZ2,MT2,MT2))*(1.0*(-(Power(EL,2)*(9 - 4*Cos(2*w) + 4*Cos(4*w))*Power(Csc(w),4)*Power(Sin(theta/2.),2))))/(1.0*(96*Power(Pi,2)));
	ComplexType PVTerm65 = (B0i(bb1,MZ2,MW2,MW2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))*Power(Cot(w),4)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm66 = (B0i(dbb0,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*Power(MB2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm67 = (B0i(dbb0,Mh12,Mh12,Mh12))*(1.0*(9*Power(Csc(w),2)*(16*Power(EL,2)*Power(Mh12,2) - (Cos(theta)*Power(3*EL*Mh12*vs*Cos(theta) + EL*Mh12*vs*Cos(3*theta) + 4*MW*(Mh12 - Mh22 + d2*Power(vs,2) + (Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),3)*Sin(w),2))/Power(vs,2))))/(1.0*(4096*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm68 = (B0i(dbb0,Mh12,Mh12,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Csc(w),2)*Power(-4*EL*(2*Mh12 + Mh22)*vs*Sin(2*theta) + 2*MW*(-5*Mh12 + 5*Mh22 + 12*d2*Power(vs,2) + 3*(Mh12 - Mh22)*Cos(2*theta))*Power(Sin(theta),2)*Sin(w) + 9*(Mh12 - Mh22)*MW*Power(Sin(2*theta),2)*Sin(w),2))))/(1.0*(8192*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm69 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(Cos(theta),3)*Power(Sin(theta),2)*Power(3*(Mh12 - Mh22)*MW*Power(Cos(theta),3) + 2*EL*(Mh12 + 2*Mh22)*vs*Csc(w)*Sin(theta) + MW*Cos(theta)*(Mh12 - Mh22 + 6*d2*Power(vs,2) - 9*(Mh12 - Mh22)*Power(Sin(theta),2)),2))))/(1.0*(1024*Power(MW,2)*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm70 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Cos(theta)*Power(Mh12 + d2*Power(vs,2) - 2*Mh22*Power(Cos(theta),2) + Mh12*Cos(2*theta),2)*Power(Sin(theta),2))))/(1.0*(256*Power(Pi,2)*Power(vs,2)));
	ComplexType PVTerm71 = (B0i(dbb0,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Power(MT2,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(16*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm72 = (B0i(dbb0,Mh12,MW2,MW2))*(1.0*(Power(EL,2)*(Power(Mh12,2) + 12*Power(MW,4) - 2*Mh12*MW2)*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sin(theta/2.),2)))/(1.0*(128*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm73 = (B0i(dbb0,Mh12,MZ2,MZ2))*(1.0*(-(Power(EL,2)*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)*(Power(Mh12,2) - 2*Mh12*MW2*Power(Sec(w),2) + 12*Power(MW,4)*Power(Sec(w),4)))))/(1.0*(256*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm74 = (B0i(dbb00,0,MB2,MB2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(12*Power(Pi,2)));
	ComplexType PVTerm75 = (B0i(dbb00,0,MT2,MT2))*(1.0*(Power(EL,2)*(-1 + Cos(theta))))/(1.0*(3*Power(Pi,2)));
	ComplexType PVTerm76 = (B0i(dbb00,0,MW2,MW2))*(1.0*(-3*Power(EL,2)*(-1 + Cos(theta))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm77 = (B0i(dbb1,Mh12,MB2,MB2))*(1.0*(3*Power(EL,2)*MB2*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm78 = (B0i(dbb1,Mh12,MT2,MT2))*(1.0*(3*Power(EL,2)*Mh12*MT2*(-1 + Power(Cos(theta),3))*Power(Csc(w),2)))/(1.0*(32*Power(MW,2)*Power(Pi,2)));
	ComplexType PVTerm79 = (B0i(dbb1,Mh12,MW2,MW2))*(1.0*(-(Power(EL,2)*Mh12*(-1 + Power(Cos(theta),3))*Power(Csc(w),2))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm80 = (B0i(dbb1,Mh12,MZ2,MZ2))*(1.0*(Power(EL,2)*Mh12*(3 + 2*Cos(theta) + Cos(2*theta))*Power(Csc(w),2)*Power(Sec(w),2)*Power(Sin(theta/2.),2)))/(1.0*(64*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24 + PVTerm25 + PVTerm26 + PVTerm27 + PVTerm28 + PVTerm29 + PVTerm30 + PVTerm31 + PVTerm32 + PVTerm33 + PVTerm34 + PVTerm35 + PVTerm36 + PVTerm37 + PVTerm38 + PVTerm39 + PVTerm40 + PVTerm41 + PVTerm42 + PVTerm43 + PVTerm44 + PVTerm45 + PVTerm46 + PVTerm47 + PVTerm48 + PVTerm49 + PVTerm50 + PVTerm51 + PVTerm52 + PVTerm53 + PVTerm54 + PVTerm55 + PVTerm56 + PVTerm57 + PVTerm58 + PVTerm59 + PVTerm60 + PVTerm61 + PVTerm62 + PVTerm63 + PVTerm64 + PVTerm65 + PVTerm66 + PVTerm67 + PVTerm68 + PVTerm69 + PVTerm70 + PVTerm71 + PVTerm72 + PVTerm73 + PVTerm74 + PVTerm75 + PVTerm76 + PVTerm77 + PVTerm78 + PVTerm79 + PVTerm80;

}

ComplexType dKappahu_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(dbb0,Mh12,Mh22,Mh22))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(dbb0,Mh12,MHA2,MHA2))*(1.0*(-(Power(del2,2)*MW2*Power(Sin(w),2))))/(1.0*(64*Power(EL,2)*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2;

}

ComplexType dgZga(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (A0i(aa0,MB2))*(1.0*(-(Power(EL,3)*(1 + 2*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2))))/(1.0*(24*MW*Power(Pi,2)));
	ComplexType PVTerm2 = (A0i(aa0,MT2))*(1.0*(Power(EL,3)*(-1 + Cos(theta))*(-8 + 3*Power(Csc(w),2))*Sec(w)))/(1.0*(24*MW*Power(Pi,2)));
	ComplexType PVTerm3 = (A0i(aa0,MW2))*(1.0*(Power(EL,3)*(2 + 3*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(8*MW*Power(Pi,2)));
	ComplexType PVTerm4 = (B0i(bb0,0,MW2,MW2))*(1.0*(Power(EL,3)*MW*(-1 + Cos(theta))*Power(Csc(w),2)*Sec(w)))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm5 = (B0i(bb0,m12,MW2,MW2))*(1.0*(-(Power(EL,3)*(6*MW2 + (Mh12 + 6*MW2)*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2))))/(1.0*(16*MW*Power(Pi,2)));
	ComplexType PVTerm6 = (B0i(bb0,m22,MW2,MW2))*(1.0*(-(Power(EL,3)*MW*(-1 + Cos(theta))*Sec(w))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm7 = (B0i(bb0,m32,MB2,MB2))*(1.0*(Power(EL,3)*MB2*(1 + 2*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(24*MW*Power(Pi,2)));
	ComplexType PVTerm8 = (B0i(bb0,m32,MT2,MT2))*(1.0*(Power(EL,3)*MT2*(-1 + 4*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(12*MW*Power(Pi,2)));
	ComplexType PVTerm9 = (B0i(bb0,m32,MW2,MW2))*(1.0*(Power(EL,3)*MW*(5 + Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm10 = (B0i(bb00,0,MB2,MB2))*(1.0*(Power(EL,3)*(1 + 2*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(12*MW*Power(Pi,2)));
	ComplexType PVTerm11 = (B0i(bb00,0,MT2,MT2))*(1.0*(Power(EL,3)*(-1 + 4*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(6*MW*Power(Pi,2)));
	ComplexType PVTerm12 = (B0i(bb00,0,MW2,MW2))*(1.0*(-(Power(EL,3)*(2 + 3*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2))))/(1.0*(4*MW*Power(Pi,2)));
	ComplexType PVTerm13 = (C0i(cc0,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-(Power(EL,3)*(m12 + m22 - m32)*MB2*(-1 + Cos(theta))*(1 + 2*Cos(2*w))*Power(Csc(w),2)*Sec(w))))/(1.0*(96*MW*Power(Pi,2)));
	ComplexType PVTerm14 = (C0i(cc0,m22,m32,m12,MT2,MT2,MT2))*(1.0*(-(Power(EL,3)*(m12 + m22 - m32)*MT2*(-1 + Cos(theta))*(-1 + 4*Cos(2*w))*Power(Csc(w),2)*Sec(w))))/(1.0*(48*MW*Power(Pi,2)));
	ComplexType PVTerm15 = (C0i(cc0,m22,m32,m12,MW2,MW2,MW2))*(1.0*(Power(EL,3)*MW*(-2*m12 + 9*m22 + m32 + Mh12 + 6*MW2 - (6*m12 - 5*m22 - 5*m32 + Mh12)*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm16 = (C0i(cc00,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-(Power(EL,3)*MB2*(1 + 2*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2))))/(1.0*(6*MW*Power(Pi,2)));
	ComplexType PVTerm17 = (C0i(cc00,m22,m32,m12,MT2,MT2,MT2))*(1.0*(Power(EL,3)*MT2*(-1 + Cos(theta))*(-8 + 3*Power(Csc(w),2))*Sec(w)))/(1.0*(6*MW*Power(Pi,2)));
	ComplexType PVTerm18 = (C0i(cc00,m22,m32,m12,MW2,MW2,MW2))*(1.0*(Power(EL,3)*(4*MW2 + (Mh12 + 6*MW2)*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(4*MW*Power(Pi,2)));
	ComplexType PVTerm19 = (C0i(cc1,m22,m32,m12,MB2,MB2,MB2))*(1.0*(Power(EL,3)*m22*MB2*(1 + 2*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(12*MW*Power(Pi,2)));
	ComplexType PVTerm20 = (C0i(cc1,m22,m32,m12,MT2,MT2,MT2))*(1.0*(Power(EL,3)*m22*MT2*(-1 + 4*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(6*MW*Power(Pi,2)));
	ComplexType PVTerm21 = (C0i(cc1,m22,m32,m12,MW2,MW2,MW2))*(1.0*(-(Power(EL,3)*MW*(-1 + Cos(theta))*(3*m12 + 7*m22 - 3*m32 + (m12 + m22 - m32)*Cos(2*w))*Power(Csc(w),2)*Sec(w))))/(1.0*(32*Power(Pi,2)));
	ComplexType PVTerm22 = (C0i(cc2,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-(Power(EL,3)*(m12 + m22 - m32)*MB2*(-1 + Cos(theta))*(1 + 2*Cos(2*w))*Power(Csc(w),2)*Sec(w))))/(1.0*(48*MW*Power(Pi,2)));
	ComplexType PVTerm23 = (C0i(cc2,m22,m32,m12,MT2,MT2,MT2))*(1.0*(-(Power(EL,3)*(m12 + m22 - m32)*MT2*(-1 + Cos(theta))*(-1 + 4*Cos(2*w))*Power(Csc(w),2)*Sec(w))))/(1.0*(24*MW*Power(Pi,2)));
	ComplexType PVTerm24 = (C0i(cc2,m22,m32,m12,MW2,MW2,MW2))*(1.0*(Power(EL,3)*MW*(4*m12 + m22 - m32 + m12*Cos(2*w))*Power(Csc(w),2)*Sec(w)*Power(Sin(theta/2.),2)))/(1.0*(8*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17 + PVTerm18 + PVTerm19 + PVTerm20 + PVTerm21 + PVTerm22 + PVTerm23 + PVTerm24;

}

ComplexType dgZga_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	return FreeTerm;

}

ComplexType dghg(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(bb0,m32,MB2,MB2))*(1.0*(Alfashgg*EL*MB2*(-1 + Cos(theta))*Csc(w)))/(1.0*(MW*Pi));
	ComplexType PVTerm2 = (B0i(bb0,m32,MT2,MT2))*(1.0*(Alfashgg*EL*MT2*(-1 + Cos(theta))*Csc(w)))/(1.0*(MW*Pi));
	ComplexType PVTerm3 = (C0i(cc0,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-(Alfashgg*EL*(m12 + m22 - m32)*MB2*Csc(w)*Power(Sin(theta/2.),2))))/(1.0*(MW*Pi));
	ComplexType PVTerm4 = (C0i(cc0,m22,m32,m12,MT2,MT2,MT2))*(1.0*(-(Alfashgg*EL*(m12 + m22 - m32)*MT2*Csc(w)*Power(Sin(theta/2.),2))))/(1.0*(MW*Pi));
	ComplexType PVTerm5 = (C0i(cc00,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-4*Alfashgg*EL*MB2*(-1 + Cos(theta))*Csc(w)))/(1.0*(MW*Pi));
	ComplexType PVTerm6 = (C0i(cc00,m22,m32,m12,MT2,MT2,MT2))*(1.0*(-4*Alfashgg*EL*MT2*(-1 + Cos(theta))*Csc(w)))/(1.0*(MW*Pi));
	ComplexType PVTerm7 = (C0i(cc1,m22,m32,m12,MB2,MB2,MB2))*(1.0*(2*Alfashgg*EL*m22*MB2*(-1 + Cos(theta))*Csc(w)))/(1.0*(MW*Pi));
	ComplexType PVTerm8 = (C0i(cc1,m22,m32,m12,MT2,MT2,MT2))*(1.0*(2*Alfashgg*EL*m22*MT2*(-1 + Cos(theta))*Csc(w)))/(1.0*(MW*Pi));
	ComplexType PVTerm9 = (C0i(cc2,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-2*Alfashgg*EL*(m12 + m22 - m32)*MB2*Csc(w)*Power(Sin(theta/2.),2)))/(1.0*(MW*Pi));
	ComplexType PVTerm10 = (C0i(cc2,m22,m32,m12,MT2,MT2,MT2))*(1.0*(-2*Alfashgg*EL*(m12 + m22 - m32)*MT2*Csc(w)*Power(Sin(theta/2.),2)))/(1.0*(MW*Pi));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10;

}

ComplexType dghg_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	return FreeTerm;

}

ComplexType dghga(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	ComplexType PVTerm1 = (B0i(bb0,m12,MW2,MW2))*(1.0*(Power(EL,3)*(Mh12 + 6*MW2)*Csc(w)*Power(Sin(theta/2.),2)))/(1.0*(8*MW*Power(Pi,2)));
	ComplexType PVTerm2 = (B0i(bb0,m22,MW2,MW2))*(1.0*(-(Power(EL,3)*MW*(-1 + Cos(theta))*Csc(w))))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm3 = (B0i(bb0,m32,MB2,MB2))*(1.0*(Power(EL,3)*MB2*(-1 + Cos(theta))*Csc(w)))/(1.0*(12*MW*Power(Pi,2)));
	ComplexType PVTerm4 = (B0i(bb0,m32,MT2,MT2))*(1.0*(Power(EL,3)*MT2*(-1 + Cos(theta))*Csc(w)))/(1.0*(3*MW*Power(Pi,2)));
	ComplexType PVTerm5 = (B0i(bb0,m32,MW2,MW2))*(1.0*(Power(EL,3)*MW*(-1 + Cos(theta))*Csc(w)))/(1.0*(16*Power(Pi,2)));
	ComplexType PVTerm6 = (C0i(cc0,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-(Power(EL,3)*(m12 + m22 - m32)*MB2*Csc(w)*Power(Sin(theta/2.),2))))/(1.0*(12*MW*Power(Pi,2)));
	ComplexType PVTerm7 = (C0i(cc0,m22,m32,m12,MT2,MT2,MT2))*(1.0*(-(Power(EL,3)*(m12 + m22 - m32)*MT2*Csc(w)*Power(Sin(theta/2.),2))))/(1.0*(3*MW*Power(Pi,2)));
	ComplexType PVTerm8 = (C0i(cc0,m22,m32,m12,MW2,MW2,MW2))*(1.0*(Power(EL,3)*(6*m12 - 5*m22 - 5*m32 + Mh12)*MW*Csc(w)*Power(Sin(theta/2.),2)))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm9 = (C0i(cc00,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-(Power(EL,3)*MB2*(-1 + Cos(theta))*Csc(w))))/(1.0*(3*MW*Power(Pi,2)));
	ComplexType PVTerm10 = (C0i(cc00,m22,m32,m12,MT2,MT2,MT2))*(1.0*(-4*Power(EL,3)*MT2*(-1 + Cos(theta))*Csc(w)))/(1.0*(3*MW*Power(Pi,2)));
	ComplexType PVTerm11 = (C0i(cc00,m22,m32,m12,MW2,MW2,MW2))*(1.0*(-(Power(EL,3)*(Mh12 + 6*MW2)*Csc(w)*Power(Sin(theta/2.),2))))/(1.0*(2*MW*Power(Pi,2)));
	ComplexType PVTerm12 = (C0i(cc1,m22,m32,m12,MB2,MB2,MB2))*(1.0*(Power(EL,3)*m22*MB2*(-1 + Cos(theta))*Csc(w)))/(1.0*(6*MW*Power(Pi,2)));
	ComplexType PVTerm13 = (C0i(cc1,m22,m32,m12,MT2,MT2,MT2))*(1.0*(2*Power(EL,3)*m22*MT2*(-1 + Cos(theta))*Csc(w)))/(1.0*(3*MW*Power(Pi,2)));
	ComplexType PVTerm14 = (C0i(cc1,m22,m32,m12,MW2,MW2,MW2))*(1.0*(-(Power(EL,3)*(m12 + m22 - m32)*MW*Csc(w)*Power(Sin(theta/2.),2))))/(1.0*(8*Power(Pi,2)));
	ComplexType PVTerm15 = (C0i(cc2,m22,m32,m12,MB2,MB2,MB2))*(1.0*(-(Power(EL,3)*(m12 + m22 - m32)*MB2*Csc(w)*Power(Sin(theta/2.),2))))/(1.0*(6*MW*Power(Pi,2)));
	ComplexType PVTerm16 = (C0i(cc2,m22,m32,m12,MT2,MT2,MT2))*(1.0*(-2*Power(EL,3)*(m12 + m22 - m32)*MT2*Csc(w)*Power(Sin(theta/2.),2)))/(1.0*(3*MW*Power(Pi,2)));
	ComplexType PVTerm17 = (C0i(cc2,m22,m32,m12,MW2,MW2,MW2))*(1.0*(Power(EL,3)*m12*MW*(-1 + Cos(theta))*Csc(w)))/(1.0*(8*Power(Pi,2)));
	return FreeTerm + PVTerm1 + PVTerm2 + PVTerm3 + PVTerm4 + PVTerm5 + PVTerm6 + PVTerm7 + PVTerm8 + PVTerm9 + PVTerm10 + PVTerm11 + PVTerm12 + PVTerm13 + PVTerm14 + PVTerm15 + PVTerm16 + PVTerm17;

}

ComplexType dghga_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32)
{

	ComplexType FreeTerm = 0;
	return FreeTerm;

}

} //end namespace SMCS
