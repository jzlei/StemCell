#ifndef BCTOOL_H
#define BCTOOL_H

//Definition of functions
#define _SQUARE(x) ((x)*(x))
#define _MAX(x,y) ((x)>=(y) ? (x):(y))
#define _MIN(x,y) ((x)<=(y) ? (x):(y))
#define _Heaviside(x) ((x)>(0) ? (1.0):(0.0))

//Definition of constants
#define ComLength 100
#define StrLength 1024
#define UNITMAX 2147483647  
#define MAXCELL 10000
#define NUMREACT 0
#define NUMVAR 5
#define NUMRAND 2
#define PI 3.14159
#define POSITIVE 0

struct IMD{
	char mdcrd[ComLength];	// trajectory file.
	char cellpar[ComLength]; // cell parameter file.
    int N0;                 // Number of cells at the initial state.
	int seed;				// Seed of the random numbers.
	double dt;              // Dynamical parameters.
    double T0;              // The time of relaxation from the initial state. We need the run the program for t < T0 so that the system reach a steady state, and than change the system parameter to study the system dynamics
	double T1;              // Dynamical parameters, total time of simulation.
	int ntpr;               // Output the result, if ntpr = 0, then no output.
	int ntpx;                // Output the path result for every ntpx steps. if ntpx = 0, then no result is outputed. 
};					    // MD parameters;

struct DCell{
    int type;
                    // type = 0: unchanged(at the resting phase of during the proliferating phase);
                    // type = 1: mitosis, generate two daughter cells;
                    // type = 2: apoptosis during the proliferating phase, removed from the proliferating phase;
                    // type = 3: terminal differentiated, removed from the resting phase
                    // type = 4: enter the proliferating state
    double X1[NUMVAR];      // State of the daughter cell 1.
    double X2[NUMVAR];      // State of the doughter cell 2.
    bool ProfQ;             // States of the cell at the proliferating phase
    double age;
};

#endif /* defined (___main__) */


