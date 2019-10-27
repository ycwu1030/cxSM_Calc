#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "ModelParameters.h"
#include "clooptools.h"
#include "CouplingFunctionSM.h"
#include "CouplingFunctionSMCS.h"
#include "CouplingFunctionSTU.h"
#include "HSS_STU_ExpData.h"

using namespace std;


int main(int argc, char* argv[])
{
    ifstream input(argv[1]);
    ofstream output(argv[2]);

    int id;
    int NTC, NTN;
    const int NMAX = 10;
    double LvhTC[NMAX], LvsTC[NMAX], HvhTC[NMAX], HvsTC[NMAX], TC[NMAX], TypeTC[NMAX];
    double LvhTN[NMAX], LvsTN[NMAX], HvhTN[NMAX], HvsTN[NMAX], TN[NMAX], TypeTN[NMAX], alpha[NMAX], betaHn[NMAX];
    RealType vh, vs, mHH, mHA, theta, a1;
    RealType mu2, lam, b1, b2, d2, del2;
    string stc, stn;


    ComplexType dkappau, dkappad, dkappac, dkappas, dkappat, dkappab;
    ComplexType dkappaele, dkappatau, dkappamuon;
    ComplexType dkappaW, dkappaZ, dkappag, dkappaga, dkappaZga;

    ltini();
    ComplexType SMkappau = SM::dKappahu(Mh12,MU2,MU2);
    ComplexType SMkappad = SM::dKappahd(Mh12,MD2,MD2);
    ComplexType SMkappac = SM::dKappahc(Mh12,MC2,MC2);
    ComplexType SMkappas = SM::dKappahs(Mh12,MS2,MS2);
    ComplexType SMkappat = SM::dKappaht(Mh12,MT2,pow(Mh1+MT+10.0,2));
    ComplexType SMkappab = SM::dKappahb(Mh12,MB2,MB2);
    ComplexType SMkappaele = SM::dKappahele(Mh12,ME2,ME2);
    ComplexType SMkappamu = SM::dKappahmuon(Mh12,MM2,MM2);
    ComplexType SMkappatau = SM::dKappahtau(Mh12,ML2,ML2);
    ComplexType SMkappaW = SM::dKappahW(Mh12,30*30,30*30);
    ComplexType SMkappaZ = SM::dKappahZ(Mh12,30*30,30*30);
    ComplexType SMkappag = SM::dghg(Mh12,0.01,0.01);
    ComplexType SMkappaga = SM::dghga(Mh12,0.01,0.01);
    ComplexType SMkappaZga = SM::dgZga(Mh12,MZ2,0.01);

    double brb=0.578;
    double brtau=0.0637;
    double brc=0.0268;
    double brgg=0.0856;
    double brgaga=0.0023;
    double brW=0.216;
    double brZ=0.0267;  
    double brmu=0.000218;
    double brOther = 1.0-brb-brtau-brc-brmu-brgg-brgaga-brW-brZ;

    double kappau, kappad, kappac, kappas, kappat, kappab, kappaele, kappatau, kappamu, kappaW, kappaZ, kappag, kappaga, kappaZga;
    double muggFb, muggFc, muggFta, muggFmu, muggFW, muggFZ, muggFga;
    double muZhb, muZhc, muZhg, muZhW, muZhta, muZhZ, muZhga, muZhmu;
    double kappaWidth;
    double S,T,U;
    bool goodTN;
    int gooditn;
    // KAPPAS kappainput;
    bool goodFitting;
    bool goodLHC8;
    bool goodLHC13;
    bool goodHLLHC300;
    bool goodHLLHC3000;
    bool goodCEPC;
    bool goodILC;
    bool goodFCC;
    int DOF;
    double chi2mu,chi2STU;
    while(input)
    {
        goodTN = false;
        input >> id >> vh >> vs >> mHH >> mHA >> theta >> a1;
        input >> mu2 >> lam >> b1 >> b2 >> d2 >> del2 >> stc >> NTC;
        for (int itc = 0; itc < NTC; ++itc)
        {
            input >> TypeTC[itc] >> LvhTC[itc] >> LvsTC[itc] >> HvhTC[itc] >> HvsTC[itc] >> TC[itc];
        }
        input >> stn >> NTN;
        for (int itn = 0; itn < NTN; ++itn)
        {
            input >> TypeTN[itn] >> LvhTN[itn] >> LvsTN[itn] >> HvhTN[itn] >> HvsTN[itn] >> TN[itn];
            input >> TypeTC[itn] >> LvhTC[itn] >> LvsTC[itn] >> HvhTC[itn] >> HvsTC[itn] >> TC[itn];
            input >> alpha[itn] >> betaHn[itn]; 
            goodTN = (goodTN || (TypeTN[itn] == 1&&betaHn[itn]>0&&alpha[itn]>0));
        }
        input.ignore(999,'\n');

        if (goodTN && mHH > sqrt(Mh12)/2.0) {
            d2 = 2*(Mh12*pow(sin(theta),2)+pow(mHH*cos(theta),2)+sqrt(2)*a1/vs)/pow(vs,2);
            dkappau = SMCS::dKappahu(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,MU2,MU2);
            dkappad = SMCS::dKappahd(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,MD2,MD2);
            dkappac = SMCS::dKappahc(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,MC2,MC2);
            dkappas = SMCS::dKappahs(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,MS2,MS2);
            dkappat = SMCS::dKappaht(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,MT2,pow(Mh1+MT+10,2));
            dkappab = SMCS::dKappahb(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,MB2,MB2);
            dkappaele = SMCS::dKappahele(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,ME2,ME2);
            dkappamuon = SMCS::dKappahmuon(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,MM2,MM2);
            dkappatau = SMCS::dKappahtau(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,ML2,ML2);
            dkappaW = SMCS::dKappahW(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,30*30,30*30);
            dkappaZ = SMCS::dKappahZ(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,30*30,30*30);
            dkappag = SMCS::dghg(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,0.01,0.01);
            dkappaga = SMCS::dghga(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,0.01,0.01);
            dkappaZga = SMCS::dgZga(mHH*mHH,mHA*mHA,theta,vs,d2,Mh12,MZ2,0.01);
            

            kappau = Re(dkappau/SMkappau)+1.0;
            kappad = Re(dkappad/SMkappad)+1.0;
            kappac = Re(dkappac/SMkappac)+1.0;
            kappas = Re(dkappas/SMkappas)+1.0;
            kappat = Re(dkappat/SMkappat)+1.0;
            kappab = Re(dkappab/SMkappab)+1.0;
            kappatau = Re(dkappatau/SMkappatau)+1.0;
            kappamu = Re(dkappamuon/SMkappamu)+1.0;
            kappaele = Re(dkappaele/SMkappaele)+1.0;
            kappaW = Re(dkappaW/SMkappaW)+1.0;
            kappaZ = Re(dkappaZ/SMkappaZ)+1.0;
            kappag = Re(dkappag/SMkappag)+1.0;
            kappaga = Re(dkappaga/SMkappaga)+1.0;
            kappaZga = Re(dkappaZga/SMkappaZga)+1.0;

            KAPPAS kappainput = {kappau,kappad,kappac,kappas,kappat,kappab,kappaele,kappamu,kappatau,kappaW,kappaZ,kappag,kappaga,kappaZga};
            kappaWidth = pow(kappab,2)*brb+pow(kappac,2)*brc+pow(kappatau,2)*brtau+pow(kappamu,2)*brmu+pow(kappaW,2)*brW+pow(kappaZ,2)*brZ+pow(kappag,2)*brgg+pow(kappaga,2)*brgaga+brOther;
            muZhb = pow(kappaZ*kappab,2)/kappaWidth;
            muZhc = pow(kappaZ*kappac,2)/kappaWidth;
            muZhg = pow(kappaZ*kappag,2)/kappaWidth;
            muZhW = pow(kappaZ*kappaW,2)/kappaWidth;
            muZhta = pow(kappaZ*kappatau,2)/kappaWidth;
            muZhZ = pow(kappaZ*kappaZ,2)/kappaWidth;
            muZhga = pow(kappaZ*kappaga,2)/kappaWidth;
            muZhmu = pow(kappaZ*kappamu,2)/kappaWidth;

            // double chisq = pow((muZhb-1)/0.0027,2)+pow((muZhc-1)/0.033,2)+pow((muZhg-1)/0.013,2)+pow((muZhW-1)/0.01,2)+pow((muZhta-1)/0.008,2)+pow((muZhZ-1)/0.051,2)+pow((muZhga-1)/0.068,2)+pow((muZhmu-1)/0.17,2);

            S = Re(STU::SinSMCS(mHH*mHH,mHA*mHA,theta,vs,d2));
            T = Re(STU::TinSMCS(mHH*mHH,mHA*mHA,theta,vs,d2));
            U = Re(STU::UinSMCS(mHH*mHH,mHA*mHA,theta,vs,d2));

            HiggsSignalStrengthSTU_Test(muLHC8,STULHC,kappainput,S,T,U,chi2mu,chi2STU,DOF,goodLHC8);
            HiggsSignalStrengthSTU_Test(muATLAS13+muCMS13,STULHC,kappainput,S,T,U,chi2mu,chi2STU,DOF,goodLHC13);
            HiggsSignalStrengthSTU_Test(muHLLHC300,STULHC,kappainput,S,T,U,chi2mu,chi2STU,DOF,goodHLLHC300);
            HiggsSignalStrengthSTU_Test(muHLLHC3000,STULHC,kappainput,S,T,U,chi2mu,chi2STU,DOF,goodHLLHC3000);
            HiggsSignalStrengthSTU_Test(muCEPC,STUCEPC,kappainput,S,T,U,chi2mu,chi2STU,DOF,goodCEPC);
            HiggsSignalStrengthSTU_Test(muILC,STUILC,kappainput,S,T,U,chi2mu,chi2STU,DOF,goodILC);
            HiggsSignalStrengthSTU_Test(muFCC,STUFCC,kappainput,S,T,U,chi2mu,chi2STU,DOF,goodFCC);
            for (int itn = 0; itn < NTN; ++itn)
            {
                if (TypeTN[itn]==1)
                {
                    gooditn = itn;
                    break;
                }
            }
            output << id << "  " << vh << "  " << vs << "  " << mHH << "  " << mHA << "  " << theta << "  " << a1 << "  " ;
            output << LvhTN[gooditn] << "  " << LvsTN[gooditn] << "  " << HvhTN[gooditn] << "  " << HvsTN[gooditn] << "  " << TN[gooditn] << "  " << LvhTN[gooditn]/TN[gooditn] << "  ";
            output << LvhTC[gooditn] << "  " << LvsTC[gooditn] << "  " << HvhTC[gooditn] << "  " << HvsTN[gooditn] << "  " << TC[gooditn] << "  " << LvhTC[gooditn]/TC[gooditn] << "  ";
            output << alpha[gooditn] << "  " << betaHn[gooditn] << "  ";
            output << kappab << "  " << kappac << "  " << kappatau << "  " << kappamu << "  ";
            output << kappaW << "  " << kappaZ << "  " << kappag << "  " << kappaga << "  " << kappaWidth << "  ";
            // output << muggFb << "  " << muggFc << "  " << muggFta << "  " << muggFmu << "  " << muggFW << "  " << muggFZ << "  " << muggFga << "  ";
            output << muZhb << "  " << muZhc << "  " << muZhg << "  " << muZhW << "  " << muZhta << "  " << muZhZ << "  " << muZhga << "  " << muZhmu << "  ";
            output << S << "  " << T << "  " << U << "  ";
            // output << chi2mu << "  " << chi2STU << "  " << DOF << "  " << goodFitting << endl;
            output << goodLHC8 << "  " << goodLHC13 << "  " << goodHLLHC300 << "  " << goodHLLHC3000 << "  " << goodCEPC << "  " << goodILC << "  " << goodFCC<< endl;
        }
        char first;
        input.get(first);
        input.putback(first);
    }
    ltexi();
    return 0;
}
