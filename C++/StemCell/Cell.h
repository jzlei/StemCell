#ifndef CCELL_H
#define CCELL_H

class CCell{
private:

	double _A[NUMVAR], _B1[NUMVAR][NUMRAND], _B2[NUMVAR][NUMRAND], _Y[NUMVAR], _Eta[NUMRAND];		// The variables to perform the stochastic simulation.
	double _Prop[NUMREACT];
  
	virtual void ReadPar(char fpar[]);					
	// Begin Gillespie algorithm
	virtual void Propensities(double *a);
	virtual void Update(int mu,double tau);
	// End of Gillespie algorithm	
	
	// Begin Stochastic simulation
	virtual void GetWienerProcess(double eta[]);
	virtual void GetRandForces(double b[][NUMRAND], double x[]);
	virtual void GetRandForces(double b[][NUMRAND], double x[], int k);
	virtual void GetTrends(double a[], double x[]);
	// End Stochastic simulation				

				
public:
	double *_X;
	int _NumReact;			// Number of reaction
	int _NumVar;			// Number of variables
	int _NumRand;			// Number of random variables
	double _t;      		// Life time after simulation
	int _cellid;
    FILE *_fp;     			// The file to record the simulation results.

	
	CCell();
	virtual ~CCell();
	
	virtual bool Initialized(int k, char fpar[]);

	virtual void GeneRead(char fpar[]);
	
	void Gillespie_OneStep();
	void Stochastic_OneStep();
	void ODE_OneStep();
	void PreRun(int type);
	void PostRun(int type);
	void PreOutPut(char fn[]);
	void Pre_OneStep();
	void Post_OneStep();
	virtual void Run(); 
	virtual void OutPut(double t);
	virtual void OutPut();
	virtual void ParOutPut(FILE *fp);
	void PostOutPut();
    
  friend class CSystem;
};

#endif