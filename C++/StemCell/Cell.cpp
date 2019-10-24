#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include <iostream>
using namespace std;

#include "BCTool.h"
#include "Random.h"
#include "Cell.h"

extern struct IMD _MD;
extern CRandom Rand;

CCell::CCell()
{
}

CCell::~CCell()
{
	delete _X;
}

void CCell::Gillespie_OneStep()
{
	double r1, r2, tau, s, s0;
	int mu;
	
	
	Propensities(_Prop);
	// Generate random pairs (\tau, \mu)
	r1=Rand();
	r2=Rand();
	tau = (1/_Prop[0])*log(1/r1);
	s=0.0;
	s0 = r2*_Prop[0];
	for (mu = 1; mu<=_NumReact; mu++)
	{
		if(s0 > s && s0 <= s + _Prop[mu])
		{
			break;
		}
		else
		{
			s=s+_Prop[mu];
		}
	}
	
	if(mu>_NumReact) { tau = 0;}	
	Update(mu,tau);   // A changed in this function is make in May 7, 2012 by JZ
}

void CCell::ODE_OneStep()
{
	int j;
	
	GetTrends(_A,_X);
	
	for(j=0;j<_NumVar;j++)
	{
		_X[j] += _A[j]*_MD.dt;
//		_X[j] = _MAX(_X[j], 0.0);
	}	
	_t += _MD.dt;
}


void CCell::Stochastic_OneStep()
{
	double f;
	int j,k;
	int j1, k1;
	
	
	GetWienerProcess(_Eta);
	GetRandForces(_B1,_X);
	GetTrends(_A,_X);
	

	// y_{t,k}^j = x_t^j + b_k^j *(eta_k^2 - dt)
	// B2_k^j = b_k^j(y_{t,k})
	for(k=0;k<_NumRand;k++){
		for(j=0;j<_NumVar;j++)
			_Y[j] = _X[j] + _B1[j][k]*((_Eta[k])*(_Eta[k]) - _MD.dt);
		GetRandForces(_B2, _Y, k);
	}
	

	for(j=0;j<_NumVar;j++){
		f = 0.0;
		for(k=0;k<_NumRand;k++){
			f += _B1[j][k]*(_Eta[k]) + 0.5*(_B2[j][k] - _B1[j][k]);}
		_X[j] += _A[j]*_MD.dt + f;
//		if (POSITIVE==1) _X[j] = _MAX(_X[j], 0.0);
	}	
	_t += _MD.dt;
}

void CCell::PreRun(int type)
{
}

void CCell::PostRun(int type)
{
}

void CCell::Run()
{
}

void CCell::GetWienerProcess(double eta[])
{
	int k;
	for(k=0;k<_NumRand;k++)
			eta[k]=Rand.WienerGen()*sqrt(_MD.dt);
}

void CCell::GetRandForces(double b[][NUMRAND], double x[])
{
	int j,k0;
	double B0[NUMVAR][NUMRAND];
	for(k0=0; k0<_NumRand; k0++)
	{
		GetRandForces(B0, x, k0);
		for(j=0; j<_NumVar; j++)
			b[j][k0]=B0[j][k0];
	}
}

void CCell::GetRandForces(double b[][NUMRAND], double x[], int k)
{
}


void CCell::GetTrends(double a[], double x[])
{
}


void CCell::Update(int mu,double tau)
{
}

bool CCell::Initialized(int k, char fpar[])
{
    return true;
}

void CCell::ReadPar(char fpar[])
{
}

void CCell::PreOutPut(char fn[])
{
     char opf[ComLength];
     if(_MD.ntpx == 0) { return;}
     sprintf(opf,"%s-%d.dat",fn,_cellid);
     if((_fp = fopen(opf,"w"))==NULL)
	 {
		 cout<<"Cannot open the output file!"<<endl;
		 exit(0);
	 }
}

void CCell::Propensities(double a[])
{
}

void CCell::OutPut(double t)
{
}

void CCell::OutPut()
{
}

void CCell::PostOutPut()
{
     if(_MD.ntpx>0) 
     {
      fclose(_fp);
     } 
}

void CCell::Pre_OneStep()
{
}

void CCell::Post_OneStep()
{
}

void CCell::GeneRead(char fpar[])
{
}

void CCell::ParOutPut(FILE *fp)
{
}



