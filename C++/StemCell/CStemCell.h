#ifndef CSTEMCELL_H
#define CSTEMCELL_H

#include "Cell.h"

struct PARStemCell{
    double beta0, kappa0, a1, a2, a3, a4, b1, eta;
    double mu0, tau;
    double theta0;
    double alpha1,alpha2;
};   // System parameters


class CStemCell:public CCell{

private:

    PARStemCell _par;
    bool _ProfQ;                       // ProfQ = 1, the cell is at the state of proliferating phase, ProfQ = 0, the cell is not at the proliferating phase (resting phase)
    double _age;                           // When the cell is at the proliferating phase, the timing _age starting from entry the proliferating rate. The cell will undergo mitosis when _age > _par.tau 

    double GetnextEpi(int i, double t, double x1, double x2);      // Get the epigenetic state of the daughter cell given the value y after DNA duplication.
    double fdeathrate();          // Get the death rate of cells.
    double fbeta(long N0, double x1, double x2);                      // The proliferation rate
    double fkappa(double u);                      // The differentiation rate

	void ReadPar(char fpar[]);
    void SetDefaultPar();
    int GetCellType(double x);     // Get the type of a cell.
    
	// Begin Stochastic simulation
	void GetRandForces(double b[][NUMRAND], double x[],int k);
	void GetTrends(double a[], double x[]);
	// End Stochastic simulation	
	
	// Begin Gillespie algorithm
	void Propensities(double a[]);
	void Update(int mu);
	// End of Gillespie algorithm
    
    void OutPutParameters();
    
	
public:


	CStemCell();
	~CStemCell();
	
	bool Initialized(int k, char fpar[]);
    DCell CellFateDecision(long N0, double dt, double t);
    
    void Pre_OneStep();
    
    
	friend class CSystem;
};

#endif
