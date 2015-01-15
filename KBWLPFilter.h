/*
 *
 *  KTransFunc
 *
 *  Copyright (c) 2000 Ryo Kikuuwe
 *
 */

#ifndef Headder_KBWFilter
#define Headder_KBWFilter

#include "KTransFunc.h"

class KBWLPFilter
{
	public:
	static void Make(KTransFunc* pTF, int Order,double aDeltaT ,double aHerz)
	{
	 //	pTF->Construct(((Order==1)?0:Order) , Order);
	 	pTF->Construct(Order , Order);
		if(aHerz<=0)
		{
			pTF->mB[0] = 1. ;
			for(int p=1;p<=pTF->N;p++) pTF->mB[p]=0.;
			for(int p=1;p<=pTF->M;p++) pTF->mA[p]=0.;
			pTF->Reset();
			return;
		}
		//if(0.5 < aHerz*aDeltaT) logMsg("[Wrn] Cutoff Freq Too Large.\n",0,0,0,0,0,0);
		double tW=2./ aDeltaT * tan( aDeltaT * PI * aHerz );
		double tG=tW * aDeltaT /2;
		if(Order==1)
		{
		//  Integral Approximation
		//	pTF->mB[0] = 1.- exp(-2. * PI * aHerz * aDeltaT )  ;
		//	pTF->mA[1] =   - exp(-2. * PI * aHerz * aDeltaT ) ;
			double tW=2./ aDeltaT * tan( aDeltaT * PI * aHerz );
			double tDenomi = 1./(tW + 2./aDeltaT) ;
			pTF->mB[0] = tDenomi * tW;
			pTF->mB[1] = pTF->mB[0];
			pTF->mA[1] = tDenomi * (tW-2./aDeltaT);
		}
		if(Order==2)
		{
			// Bilinear Transform
			double tDenomi = 1./(1. + sqrt(2.)*tG + tG*tG );
			pTF->mB[0] = tDenomi * tG*tG;
			pTF->mB[1] = 2.* pTF->mB[0];
			pTF->mB[2] = pTF->mB[0];
			pTF->mA[1] = tDenomi * 2.*(tG*tG-1.);
			pTF->mA[2] = tDenomi * (tG*tG - sqrt(2.)*tG +1. );
		}
		if(Order==3)
		{
			double tDenomi = 1./(1. + tG)/(1. + tG + tG * tG ) ;
			pTF->mB[0]=tDenomi * ( tG*tG*tG )  ;
			pTF->mB[1]= 3.* pTF->mB[0] ;
			pTF->mB[2]= 3.* pTF->mB[0] ;
			pTF->mB[3]= pTF->mB[0] ;
			pTF->mA[1]=tDenomi * ( -3. - 2.*tG + 2.*tG*tG+ 3.*tG*tG*tG ) ;
			pTF->mA[2]=tDenomi * ( 1. + tG )*(3. - 5.*tG + 3.*tG*tG )  ;
			pTF->mA[3]=tDenomi * ( -1. + 2.*tG - 2.*tG*tG + tG*tG*tG ) ;
		}
		if(Order==4)
		{
			double tDenomi = 1./(1. + 2.*cos(PI/8)*tG + tG*tG )/(1. + 2.*cos(3.*PI/8)*tG + tG*tG) ;
			pTF->mB[0]= tDenomi * tG*tG*tG*tG ;
			pTF->mB[1]= pTF->mB[0]*4. ;
			pTF->mB[2]= pTF->mB[0]*6. ;
			pTF->mB[3]= pTF->mB[0]*4. ;
			pTF->mB[4]= pTF->mB[0] ;
			pTF->mA[1]= tDenomi* 4.*(-1. + tG*tG)*(1. + (cos(PI/8) +cos(3.*PI/8)) *tG + tG*tG) ;
			pTF->mA[2]= tDenomi* 2.*( 3. - (2.+sqrt(2.))*tG*tG + 3. *tG*tG*tG*tG) ;
			pTF->mA[3]= tDenomi* 4.*(-1. + tG*tG)*(1. - (cos(PI/8)+cos(3*PI/8)) * tG  + tG*tG)   ;
			pTF->mA[4]= tDenomi* (1. - 2.*cos(PI/8)*tG + tG*tG) * (1. - 2.*tG*cos(3.*PI/8) + tG*tG)  ;
	  	}
		pTF->Reset();
	}
	//----------------------------------------------------
	static void QuickMake(KTransFunc* pTF, int Order, double aDeltaT ,double aHerz)
	{
		if(aHerz<=0)
		{
			pTF->mB[0] = 1. ;
			for(int p=1;p<=pTF->N;p++) pTF->mB[p]=0.;
			for(int p=1;p<=pTF->M;p++) pTF->mA[p]=0.;
		}
		if(Order==1)
		{
		// Integral Approximation
		//	pTF->mB[0] = 1.- exp(-2. * PI * aHerz * aDeltaT )  ;
		//	pTF->mA[1] =   - exp(-2. * PI * aHerz * aDeltaT ) ;
		// H = T/(exp(2 aHerz PI*aDeltaT)-1.);
		//	double HH = 1./(2.*PI*aHerz) ;
		//	pTF->mB[0] = aDeltaT/( HH + aDeltaT ) ;
		//	pTF->mA[1] = HH     /( HH + aDeltaT ) ;
			double tW=2./ aDeltaT * tan( aDeltaT * PI * aHerz );
			double tDenomi = 1./(tW + 2./aDeltaT) ;
			pTF->mB[0] = tDenomi * tW;
			pTF->mB[1] = pTF->mB[0];
			pTF->mA[1] = tDenomi * (tW-2./aDeltaT);
		}
		if(Order==2)
		{
			// Bilinear Transpose
			double tW=2./ aDeltaT * tan( aDeltaT * PI * aHerz );
			double tG=tW * aDeltaT /2;
			double tDenomi = 1./(1. + sqrt(2.)*tG + tG*tG );
			pTF->mB[0] = tDenomi * tG*tG;
			pTF->mB[1] = 2.* pTF->mB[0];
			pTF->mB[2] = pTF->mB[0];
			pTF->mA[1] = tDenomi * 2.*(tG*tG-1.);
			pTF->mA[2] = tDenomi * (tG*tG - sqrt(2.)*tG +1. );
		}
	}
	static void QuickCopy(KTransFunc* pTF, const KTransFunc& aTF)
	{
      {for(int i=0;i<=pTF->N;i++) pTF->mB[i]=aTF.mB[i];}
      {for(int i=0;i<=pTF->M;i++) pTF->mA[i]=aTF.mA[i];}
	}
	//----------------------------------------------------
///////////////////////////////////////////////////////////////
};


#endif
