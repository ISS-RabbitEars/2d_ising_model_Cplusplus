#include "2DIsing.h"
#include <sstream>
using namespace std;

struct itr 
{
    double o,f,i;
    itr operator=(double a) {o=f=i=a; return *this;}
    friend istream& operator>>(istream &in,itr &a) {in>>a.o; in>>a.f; in>>a.i; return in;}
};

struct array2D
{
    int r,c;
    double **a;
    array2D() {r=c=0;}
    ~array2D() {for(int i=0;i<r;i++) {delete[] a[i];} delete[] a; a=NULL;}
    void init(int x,int y) {r=x;c=y; a=new double*[x]; for(int i=0;i<r;i++) {a[i]=new double[c];} for(int i=0;i<r;i++) for(int j=0;j<c;j++) a[i][j]=0.;}
};

void ReadInItr(itr&,string&);

int main(int argc, char* argv[])
{
    ising2d sim;
    itr T,B; T=0.; B=0.;
    string pfn=argv[1], tfn=argv[2], bfn=argv[3],rootfnT,rootfnB,tag,up,down;
    array2D uU,uM,uC,uX,dU,dM,dC,dX;
    int x,y;

    rootfnT="T=";
    rootfnB="B=";
    tag=".dat";
    up="_ramp_up";
    down="_ramp_down";
    ReadInItr(T,tfn);
    ReadInItr(B,bfn);
    x=int((T.f-T.o)/T.i)+1;
    y=int((B.f-B.o)/B.i)+1;
    uU.init(x,y);uM.init(x,y);uC.init(x,y);uX.init(x,y);dU.init(x,y);dM.init(x,y);dC.init(x,y);dX.init(x,y);
    sim.init(pfn.c_str());
    int i=0;
    for(double t=T.o;t<=T.f;t+=T.i)
    {
        int j=0;
        for(double b=B.o;b<=B.f;b+=B.i)
        {
            sim.SMS(t,b);
            uU.a[i][j]=sim.avgu(); uM.a[i][j]=sim.avgM(); uC.a[i][j]=sim.avgc(); uX.a[i][j]=sim.avgX();
            j++;
        }
        j=0;
        for(double b=B.f;b>=B.o;b-=B.i)
        {
            sim.SMS(t,b);
            dU.a[i][j]=sim.avgu(); dM.a[i][j]=sim.avgM(); dC.a[i][j]=sim.avgc(); dX.a[i][j]=sim.avgX();
            j++;
        }
        i++;
    }
    
    string fu="_u(T)",fM="_M(T)",fC="_C(T)",fX="_X(T)";
    ofstream foutuU,foutuM,foutuC,foutuX,foutdU,foutdM,foutdC,foutdX;
    for(i=0;i<y;i++)
    {
        string fnuU,fnuM,fnuC,fnuX,fndU,fndM,fndC,fndX,snum;
        stringstream ssnum;
        ssnum<<B.o+(i*B.i);
        snum=ssnum.str();
        ssnum.str("");
        fnuU=rootfnB+snum+up+fu+tag;
        fndU=rootfnB+snum+down+fu+tag;
        fnuM=rootfnB+snum+up+fM+tag;
        fndM=rootfnB+snum+down+fM+tag;
        fnuC=rootfnB+snum+up+fC+tag;
        fndC=rootfnB+snum+down+fC+tag;
        fnuX=rootfnB+snum+up+fX+tag;
        fndX=rootfnB+snum+down+fX+tag;
        foutuU.open(fnuU.c_str());
        foutdU.open(fndU.c_str());
        foutuM.open(fnuM.c_str());
        foutdM.open(fndM.c_str());
        foutuC.open(fnuC.c_str());
        foutdC.open(fndC.c_str());
        foutuX.open(fnuX.c_str());
        foutdX.open(fndX.c_str());
        for(int j=0;j<x;j++)
        {
            double Ti=T.o+(j*T.i);
            foutuU<<Ti<<" "<<uU.a[j][i]<<endl;
            foutdU<<Ti<<" "<<dU.a[j][y-1-i]<<endl;
            foutuM<<Ti<<" "<<uM.a[j][i]<<endl;
            foutdM<<Ti<<" "<<dM.a[j][y-1-i]<<endl;
            foutuC<<Ti<<" "<<uC.a[j][i]<<endl;
            foutdC<<Ti<<" "<<dC.a[j][y-1-i]<<endl;
            foutuX<<Ti<<" "<<uX.a[j][i]<<endl;
            foutdX<<Ti<<" "<<dX.a[j][y-1-i]<<endl;
        }
        foutuU.close();
        foutdU.close();
        foutuM.close();
        foutdM.close();
        foutuC.close();
        foutdC.close();
        foutuX.close();
        foutdX.close();
    }
    
    return 0;
}

void ReadInItr(itr &b, string &fn)
{
    ifstream fin;
    fin.open(fn.c_str());
    fin>>b;
    fin.close();
}