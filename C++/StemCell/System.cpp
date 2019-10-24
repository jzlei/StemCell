#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include <iostream>
using namespace std;

#include "BCTool.h"
#include "Cell.h"
#include "CStemCell.h"
#include "System.h"
#include "Random.h"

extern struct IMD _MD;
extern CRandom Rand;

CSystem::CSystem()
{
	_NumCell = 0;
}

CSystem::CSystem(int N0)
{
    _Prolif=1.0;
    
    _NumPoolCell = N0;      // Number of cells in the simulation pool
    _NumCell = N0;          // Number of total cell numbers
    _MaxNumCell = MAXCELL;
    
    _cells = new CStemCell[_MaxNumCell];
}

CSystem::~CSystem()
{
     delete _cells;
}

bool CSystem::Initialized()
{
	int k;
	k=1;
	do
	{
		if((*this)(k).Initialized(k, _MD.cellpar))
		{
			k++;
		}
		else
		{
		}
	}while(k<=_MaxNumCell);

	return(true);
}

bool CSystem::SystemUpdate(double t)
{
    int k;
    int Ntemp;
    double *X1, *X2;      //  X1 stores the state x1, and X2 stores x2 for all cells
    bool *PQ;     //  store the state of _ProlifQ for all cells
    double *age;    //  store the state of age for all cells
    
    X1 = new double[2*_NumPoolCell+1];
    X2 = new double[2*_NumPoolCell+1];
    PQ = new bool[2*_NumPoolCell+1];
    age = new double[2*_NumPoolCell+1];
    Ntemp = 0;        // Number of cells after a cycle.
    
    _N0 = 0;          // Number of cells in the pool after cell fate decision.
    _N1 = 0;          // Number of cells removed from resting phase
    _N2 = 0;          // Number of cells undergoing mitosis
    _N3 = 0;          // Number of cells removed from the proliferating phase
    _N4 = 0;          // Number of cells remain unchanged at the proliferating phase
    _N5 = 0;          // Number of cells remain unchanged at the resting phase
    
    
    for (k=1; k<=_NumPoolCell; k++)
    {
        nextcell=(*this)(k).CellFateDecision(_NumCell, _MD.dt, t);
        switch(nextcell.type){
            case 3:                 // Cells with ternimal differentiation
                _N1++;
                break;
            case 0:                 // Cells remain unchanged
            case 4:                 // Enter the prolifertaing state
                Ntemp++;
                X1[Ntemp] = nextcell.X1[0];
                X2[Ntemp] = nextcell.X1[1];
                PQ[Ntemp] = nextcell.ProfQ;
                age[Ntemp] = nextcell.age;
                _N0++;
                if(nextcell.ProfQ==0)
                {
                    _N5++;
                }
                else
                {
                    _N4++;
                }
                
                break;
            case 1:                 // Mitosis
                Ntemp++;
                X1[Ntemp] = nextcell.X1[0];
                X2[Ntemp] = nextcell.X1[1];
                PQ[Ntemp] = nextcell.ProfQ;
                age[Ntemp] = nextcell.age;
                
                Ntemp++;
                X1[Ntemp] = nextcell.X2[0];
                X2[Ntemp] = nextcell.X2[1];
                PQ[Ntemp] = nextcell.ProfQ;
                age[Ntemp] = nextcell.age;
                

                _N0 = _N0+2;
                _N2 = _N2+1;
                
                break;
            case 2:                 // Apoptosis during proliferation
                _N3++;
                break;
        }
        
    }
    
    _Prolif = 1.0*Ntemp/_NumPoolCell;
    _NumCell =_NumCell * _Prolif;
    
    if(Ntemp==0)
    {
        return(0);
    }
    
    int i;
    double p0;

    p0 = 1.0*_MaxNumCell/Ntemp;
    k=0;
    for (i=1; i<=Ntemp; i++)
    {
        if(Rand()<p0 && k<_MaxNumCell)
        {
            k=k+1;
            (*this)(k)._X[0] = X1[k];
            (*this)(k)._X[1] = X2[k];
            (*this)(k)._ProfQ = PQ[k];
            (*this)(k)._age = age[k];
        }
        if (k==_MaxNumCell)
        {
            break;
        }
    }
    _NumPoolCell = k;
    
    delete X1;
    delete X2;
    delete PQ;
    delete age;
    return(1);
}

void CSystem::Run()
{
    double t;
    int step;
    int k;
    FILE *fmdsys;
    char fn[ComLength];
    sprintf(fn,"%s.dat",_MD.mdcrd);
    if((fmdsys = fopen(fn,"w"))==NULL)
    {
        cout<<"Cannot open the file mdsys."<<endl;
        exit(0);
    }
    
    step=0;
    OutPutSys(step);

    for (t=0; t<=_MD.T1; t=t+_MD.dt)
    {
        step=step+1;
        for (k=1; k<=_NumPoolCell; k++){
            (*this)(k).Pre_OneStep();
        }
        if(SystemUpdate(t))
        {
            fprintf(fmdsys,"%f %f %f %d %d %d %d %d %d\n",t,_Prolif, _NumCell, _NumPoolCell, _N0, _N1, _N2, _N3, _N4);
            if((_MD.ntpx>0) && (step%_MD.ntpx==0))
            {
                OutPutSys(step);
            }
        }
        else{
            printf("Tumor cells cleaned at day %4.2f\n",t);
            break;
        }
        
    }
    fclose(fmdsys);
}

void CSystem::RunOneStep(double t)
{
}

void CSystem::OutPutSys(int step)
{
    FILE *fp;
    int k;
    char fnc[StrLength];
    
    sprintf(fnc,"%s-%d.dat",_MD.mdcrd,step);
    if((fp = fopen(fnc,"w"))==NULL)
    {
        cout<<"Cannot open the file fmdc."<<endl;
        exit(0);
    }

	
	for(k = 1; k<= _NumPoolCell; k++)
	{
	// Edit to change you output items.

        fprintf(fp,"%d %6.3f %6.3f\n",k, (*this)(k)._X[0], (*this)(k)._X[1]);
	}
	 
    fclose(fp);
}

void CSystem::OutPutCells(double t)
{
	int k;
	for(k = 1; k<= _NumPoolCell; k++)
	{
        (*this)(k).OutPut(t);
	}
}


CStemCell& CSystem::operator()(int indx)
{
if(indx>0 && indx<=_MaxNumCell)
{
	return *(_cells+(indx-1));
}
else
{
	cout<<"Err<< CSystem () >>Dimensional error"<<endl;
	exit(0);
}
}
