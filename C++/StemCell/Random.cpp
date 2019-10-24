#include "Random.h"

#define RANDMAX 23371.0
#define RANDPI 3.1415926

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

CRandom::CRandom()
{
    int i;
    time_t t;
    double f;
    srand( (unsigned) time(&t));
    for (i = 0; i < 97; i++)
    {
        f = 1.0*rand();
        x[i] = (1.0*fmod(f,RANDMAX))/RANDMAX;
    }
    c = (1.0*fmod(rand(),RANDMAX))/RANDMAX;
}

CRandom::CRandom(unsigned seed)
{
    Initialized(seed);
}

void CRandom::Initialized(unsigned seed)
{
    int i;
    double f;
    srand(seed);
    for (i = 0; i < 97; i++)
    {
        f = 1.0*rand();
        x[i] = (1.0*fmod(f,RANDMAX))/RANDMAX;
    }
    c = (1.0*fmod(rand(),RANDMAX))/RANDMAX;
}

double CRandom::Value()
{
    double x0,U;
    int r = 0, s = 64;
    int i;
    double d = 7654321.0/16777216.0, d0 = 1677213.0/1677216.0;
    double f;
    if (x[r] >= x[s])
    {
        x0 = x[r] - x[s];
    }
    else
    {
        x0 = x[r] - x[s] + 1;
    }
    if ( c >= d)
    {
        c = c - d;
    }
    else
    {
        c = c - d + d0;
    }
    if (x0 >= c)
    {
        U = x0 - c;
    }
    else
    {
        U = x0 - c + d0;
    }
    for (i = 0; i < 96; i++)
    {
        x[i] = x[i+1];
    }
    x[96] = fmod(x0,1);
    c = fmod(c,1);
    f = fmod(U,1);
    return(f);
}

// Get the value with standard normal distribution;
double CRandom::NormalDistribution()
{
    double f1,f2,x;
    f1 = Value();
    f2 = Value();
    x = sqrt(-2*log(f1))*cos(2*RANDPI*f2);
    return(x);
}

// Get the value with normal distribution with mean mu and standard deviation sigma
double CRandom::NormalDistribution(double mu, double sigma)
{
    double x;
    x = NormalDistribution();
    x = sigma*x+mu;
    return(x);
}

// Get a random number with Binormal distribution B(n,p)
int CRandom::BinomialDistribution(int n, double p)
{
    double r;
    double a0, a1, sum;
    int k;
    r = Value();
    if (p > 1.0 - 1e-6)
    {
        k=n;
    }
    else
    {
        a0 = pow(1-p,n);
        sum = a0;
        if (sum > r){
            k=0;
        }
        else{
            for (k=1; k<=n; k++)
            {
                a1 = a0 * (1.0 * (n-k + 1.0)/(1.0 * k))*(p/(1-p));
                sum = sum + a1;
                if (sum > r)
                {
                    break;
                }
                else
                {
                    a0 = a1;
                }
            }
        }
    }
    return(k);
}

double CRandom::operator()()
{
    return(Value());
}

double CRandom::operator()(double a, double b)
{
    return(Value(a,b));
}


// Obtain a random value from the interval [a,b]
double CRandom::Value(double a, double b)
{
    return(a + (b-a)*Value());
}

double CRandom::WienerGen()
{
    double f1,f2,x;
    f1 = Value();
    f2 = Value();
    x = sqrt(-2*log(f1))*cos(2*RANDPI*f2);
    return(x);
}

double CRandom::BetaDistribution(double a, double b)
{
    double x, y;
    x = GammaDistribution(a,1);
    y = GammaDistribution(b,1);
    return(x/(x+y));
}

double CRandom::GammaDistribution(double a, double b)
{
    int k,n;
    double delta, U, V, W, xi, eta;
    double E=2.71828;
    n=floor(a);
    delta = a - n;
    xi=0;
    if(delta>1e-6)
    {
        for (k=1; k<=100; k++)
        {
            U=Value();
            V=Value();
            W=Value();
            if(U<=E/(E+delta))
            {
                xi = pow(V,1/delta);
                eta = W * pow(xi, delta - 1);
            }
            else
            {
                xi = 1 - log(V);
                eta = W * pow(E,-xi);
            }
            if (eta < pow(xi, delta-1)*pow(E, -xi))
            {
                break ;
            }
        }
    }
    for (k=1;k<=n;k++)
    {
        xi = xi - log(Value());
    }
    return(b*xi);
}
