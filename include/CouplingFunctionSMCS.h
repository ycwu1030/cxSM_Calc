#ifndef CouplingFunctionSMCS
#define CouplingFunctionSMCS
#include "clooptools.h"

namespace SMCS{
ComplexType dKappahW(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahW_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappahZ(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahZ_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappahb(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahb_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappahc(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahc_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappahd(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahd_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappahele(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahele_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappahmuon(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahmuon_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappahs(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahs_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappaht(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappaht_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappahtau(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahtau_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dKappahu(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dKappahu_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dgZga(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dgZga_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dghg(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dghg_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

ComplexType dghga(double Mh22, double MHA2, double theta, double vs, double d2, double m12, double m22, double m32);

ComplexType dghga_Z2(double Mh22, double MHA2, double d2, double del2, double m12, double m22, double m32);

} //end namespace SMCS
#endif
