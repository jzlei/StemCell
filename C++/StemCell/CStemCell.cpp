#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include <iostream>
using namespace std;

#include "Random.h"
#include "BCTool.h"
#include "CStemCell.h"

extern struct IMD _MD;
extern CRandom Rand;
extern double _ratio19;


CStemCell::CStemCell()
{
}

CStemCell::~CStemCell()
{
}

bool CStemCell::Initialized(int k, char fpar[])
{
    SetDefaultPar();
	ReadPar(fpar);
	_t = 0.0;
	_cellid=k;
	_NumVar=2;
	_X = new double[_NumVar];
    _X[0] = Rand(0,1);        // Initial state of the epigenetic state x1
    _X[1] = Rand(0,0.1);      // Initial state of the epigenetic state x2
    _ProfQ = 0;               // We assume that all cells at resting phase at the initial state.
    _age = 0;
    
    return true;
}

void CStemCell::ReadPar(char fpar[])
{
	FILE *fp;
	char str[StrLength], *pst;
	if((fp = fopen(fpar,"r"))==NULL)
	{
		cout<<"Cannot open the cell parameter input file."<<endl;
		exit(0);
	}
	rewind(fp);
	while(!feof(fp))
	{
		fgets(str,StrLength,fp);
		if(str[0]=='#'){ continue;}

        if((pst=strstr(str,"beta0="))!=NULL)
        {
            _par.beta0=atof(pst+6);
        }
        if((pst=strstr(str,"kappa0="))!=NULL)
        {
            _par.kappa0=atof(pst+7);
        }
        if((pst=strstr(str,"p_a1="))!=NULL)
        {
            _par.a1=atof(pst+5);
        }
        if((pst=strstr(str,"p_a2="))!=NULL)
        {
            _par.a2=atof(pst+5);
        }
        if((pst=strstr(str,"p_a3="))!=NULL)
        {
            _par.a3=atof(pst+5);
        }
        if((pst=strstr(str,"p_a4="))!=NULL)
        {
            _par.a4=atof(pst+5);
        }
        if((pst=strstr(str,"p_b1="))!=NULL)
        {
            _par.b1=atof(pst+5);
        }
        if((pst=strstr(str,"eta="))!=NULL)
        {
            _par.eta=atof(pst+4);
        }
        if((pst=strstr(str,"mu="))!=NULL)
        {
            _par.mu0=atof(pst+3);
        }
        if((pst=strstr(str,"tau="))!=NULL)
        {
            _par.tau=atof(pst+4);
        }
        if((pst=strstr(str,"theta0="))!=NULL)
        {
            _par.theta0=atof(pst+7);
        }
        if((pst=strstr(str,"alpha1="))!=NULL)
        {
            _par.alpha1=atof(pst+7);
        }
        if((pst=strstr(str,"alpha2="))!=NULL)
        {
            _par.alpha2=atof(pst+7);
        }
 	}
	fclose(fp);
    
//    OutPutParameters();    
}

void CStemCell::OutPutParameters()
{
    printf("beta0=%f, theta0=%f\n", _par.beta0, _par.theta0);
    printf("kappa0=%f, mu=%f, tau=%f\n",_par.kappa0, _par.mu0, _par.tau);
    printf("a1=%f, a2=%f, a3=%f, a4=%f, b1=%f, alpha1=%f, alpha2=%f\n", _par.a1, _par.a2, _par.a3, _par.a4, _par.b1, _par.alpha1, _par.alpha2);
}

void CStemCell::SetDefaultPar()
{
}

void CStemCell::GetRandForces(double b[][NUMRAND], double x[], int k)
{
}

void CStemCell::GetTrends(double a[], double x[])
{
}

DCell CStemCell::CellFateDecision(long N0, double dt, double t)
{
    int i;
    double mu;
    double beta;
    double kappa;
    double _rand;
    
    mu = fdeathrate() * dt;
    beta = fbeta(N0, _X[0], _X[1]) * dt;
    kappa = fkappa(_X[0]) * dt;
    
    DCell nextcell;
    
    nextcell.type = 0;  //Default, cells remain unchanged.
    for (i=0;i<_NumVar;i++)
    {
        nextcell.X1[i] = _X[i];
    }
    nextcell.ProfQ = _ProfQ;
    nextcell.age = _age;
    
    if (_ProfQ == 0)        // If the cell is at the resting phase
    {
        _rand = Rand();
        if (_rand < kappa)  // Terminal differentiation from the resting phase
        {
            nextcell.type = 3;
        }
        else
        {
            if (_rand < kappa + beta)   // Enter the proliferting phase
            {
                nextcell.type = 4;
                for (i=0;i<_NumVar;i++)
                {
                    nextcell.X1[i] = _X[i];
                }
                nextcell.ProfQ = 1;
                nextcell.age = 0;
            }
        }
    }
    
    if (_ProfQ == 1)    // If the cell is at the proliferating phase
    {
        _rand = Rand();

        if (_rand < mu)     // The apoptosis during the proliferating phase
        {
            nextcell.type = 2;
        }
        else
        {
            if (_age < _par.tau)
            {
                nextcell.age = _age + dt;
            }
            else  // Perform mitosis
            {
                nextcell.type = 1;
                for (i=0;i<_NumVar;i++)
                {
                    nextcell.X1[i] = GetnextEpi(i, t, _X[0], _X[1]);
                    nextcell.X2[i] = GetnextEpi(i, t, _X[0], _X[1]);
                }
                nextcell.ProfQ = 0;
                nextcell.age = 0;
            }
        }
    }
    
    return(nextcell);
}

double CStemCell::GetnextEpi(int i, double t, double x1, double x2)
{
    double phi;
    double a, b;
    double z;
    double f;
    switch(i)
    {
        case 0:
            phi = 0.08 + 1.0 * pow(_par.alpha1 * x1, 1.8)/(1+pow(_par.alpha1 * x1, 1.8));
            break;
        case 1:
            f = 1.0/(1 + exp(-(t - _MD.T0)/1000.0));
            phi = 0.08 + (1.0 + f * 0.4/(1 + pow(2.5 * x1, 6))) * pow(_par.alpha2 * x2, 1.8)/(1+pow(_par.alpha2 * x2, 1.8));
            break;
    }
    a = _par.eta  * phi;
    b = _par.eta * (1-phi);
    z = Rand.BetaDistribution(a, b);
    return(z);
}

void CStemCell::Update(int mu)
{
}

void CStemCell::Propensities(double a[])
{
}

double CStemCell::fbeta(long N0, double x1, double x2)
{
    double beta, beta0;
    double theta;
    
    beta0 = _par.beta0;
    theta = _par.theta0 * (1.0 + pow(_par.a4 * x2, 6.0)/(1.0 + pow(_par.a4 * x2, 6.0)));
    
    beta = beta0 * (1.0/(1.0 + N0/theta)) * ((_par.a1 * x1 + pow(_par.a2 * x1,6.0))/(1+pow(_par.a3 * x1,6.0)));   // Proliferation rate of cells
    return(beta);
}

double CStemCell::fkappa(double x1)
{
    double kappa;
    
    kappa = _par.kappa0 * 1.0/(1.0 + pow(_par.b1 * x1, 6.0));

    return(kappa);
}

double CStemCell::fdeathrate()
{
    double mu;
    
    mu = _par.mu0;
    
    return(mu);
}

void CStemCell::Pre_OneStep()
{
}



