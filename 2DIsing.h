#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define Nbrs 4

class ising2d
{
    int N,m,n,ns,teq,mts,*s,**nn;
    double J,kB,mo,T,B;
    double rN,mM,rmts,mU,mU2,mM2;
    
public:
     
    ising2d();
    ~ising2d();
    void init(string); //read in parameters from file and initialize arrays
    void init(int,int,double,double,double); //pass parameters and initialize arrays (m,n,J,kB,mo)
    void WriteNNList();
    void ScanMicroStates(double,double); //equilibrate system for given T and B by scanning over 2^N microstates using Metropolis Monte Carlo to minimize energy
    void ScanMicroStates(double,double,int);
    void SMS(double,double);
    double SumNN(int); //nearest neighbor spin sum for particle i 
    double SumSpinProduct(); //
    double SumSpins(); //
    double U(); //calculate the total potential energy of the system
    double u(); //potential energy per particle
    double M();
    void rescale();
    double avgM();
    double avgU();
    double avgu();
    double avgC();
    double avgc();
    void initM(double&);
    void update1M(double&,int);
    void update2M(double);
    void initU(double&);
    void update1U(double&,double);
    void update2U(double);
    double avgX();

    
    double ran3(int*); 
};

#include "2DIsing.cpp"
#undef Nbrs
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
