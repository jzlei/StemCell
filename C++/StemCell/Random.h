#ifndef RANDOM_H
#define RANDOM_H

//#define RANDMAX 23371.0
//#define RANDPI 3.1415926

class CRandom{
private:

	double x[97];
    double c;
public:

	 
	CRandom();
    CRandom(unsigned seed);
	double Value();                               // Get the value with uniform distribution;
	double NormalDistribution();      // Get the value with standard normal distribution;
	double NormalDistribution(double mu, double sigma);  // Get the value with normal distribution with mean mu and standard deviation sigma
    int BinomialDistribution(int n, double p);      // Get a random number with binormal distribution B(n,p);
    double BetaDistribution(double a, double b);    // Get a random number with BetaDistribution
    double GammaDistribution(double a, double b);   // Get a random number with GammaDistribution
	void Initialized(unsigned seed);    // Initialzation of the random number

	double operator()();                   // Get a random value;
	double Value(double a, double b);
	double operator()(double a, double b);
    double WienerGen();                    // Generate a Wiener process;
};

#endif


