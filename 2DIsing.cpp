ising2d::ising2d()
{
    N=m=n=ns=teq=mts=0;
    J=kB=mo=T=B=rN=mM=mU=mU2=rmts=mM2=0.;
}

ising2d::~ising2d()
{
    if(N!=0)
    {
        delete[] s; s=NULL;
        for(int i=0;i<N;i++) delete[] nn[i];
        delete[] nn; nn=NULL;
    }
}

void ising2d::init(string fn)
{
    int np=5;
    double p[np];
    ifstream fin;
    fin.open(fn.c_str());
    for(int i=0;i<np;i++) fin>>p[i];
    fin.close();
    init(int(p[0]),int(p[1]),p[2],p[3],p[4]);
}

void ising2d::init(int x, int y, double ee, double kb, double dm)
{
    m=x; n=y; N=m*n; J=ee; kB=kb; mo=dm;
    rN=1./double(N); 
    ns=20*int(pow(double(N),2)); teq=int(0.666*double(ns)); mts=ns-teq; rmts=1./double(mts); //pow(2.,N)
    s=new int[N];
    nn=new int*[N];
    int seed=-426156387;
    for(int i=0;i<N;i++) {(ran3(&seed)>=0.5)?s[i]=-1:s[i]=1; nn[i]=new int[Nbrs];}
    int a=m-1,b=m*(n-1),c=N-1;
    nn[0][0]=m; nn[0][1]=1; nn[0][2]=b; nn[0][3]=a;
    nn[a][0]=(2*m)-1; nn[a][1]=0; nn[a][2]=c; nn[a][3]=a-1;
    nn[b][0]=0; nn[b][1]=b+1; nn[b][2]=b-m; nn[b][3]=c;
    nn[c][0]=a; nn[c][1]=b; nn[c][2]=c-m; nn[c][3]=c-1;
    for(int i=1;i<a;i++) {int j=b+i; nn[i][0]=i+m; nn[i][1]=i+1; nn[i][2]=j; nn[i][3]=i-1; nn[j][0]=i; nn[j][1]=j+1; nn[j][2]=j-m; nn[j][3]=j-1;}
    for(int i=m;i<b;i+=m) {int j=i+m-1; nn[i][0]=i+m; nn[i][1]=i+1; nn[i][2]=i-m; nn[i][3]=j; nn[j][0]=j+m; nn[j][1]=i; nn[j][2]=j-m; nn[j][3]=j-1;}
    for(int j=1;j<n-1;j++) for(int i=1;i<a;i++) {int k=(j*m)+i; nn[k][0]=k+m;nn[k][1]=k+1;nn[k][2]=k-m;nn[k][3]=k-1;}
}

void ising2d::WriteNNList()
{
     ofstream fout;
     fout.open("NN_List.dat");
     for(int i=0;i<N;i++) {fout<<i; for(int j=0;j<Nbrs;j++) {fout<<" "<<nn[i][j];} fout<<endl;}
     fout.close();
}

void ising2d::ScanMicroStates(double t, double b)
{
    T=t; B=b;
    int c=0,cn=0,seed=-845346743;
    double rkbt=1./(kB*T);
    double rns=1./double(ns);
    for(int j=0;j<ns;j++)
    {
        int i=int(double((N-1))*ran3(&seed));
        double dE=2.*double(s[i])*((J*SumNN(i))+(B*mo));
        if(dE<=0) {s[i]=-s[i]; c++;}
        else {double p=exp(-dE*rkbt); if(p>=ran3(&seed)) {s[i]=-s[i]; c++;} else cn++;}
        cout<<c<<"  "<<cn<<"  "<<double(c)/double(cn)<<"  "<<double(c+cn)*rns<<"  "<<double(c)*rns<<endl;
    }
}

void ising2d::ScanMicroStates(double t, double b, int mod)
{
    T=t; B=b; mM=0.;
    int seed=-845346743;
    double rkbt=1./(kB*T);
    for(int j=0;j<ns;j++)
    {
        int i=int(double((N-1))*ran3(&seed));
        double dE=2.*double(s[i])*((J*SumNN(i))+(B*mo));
        if(dE<=0) {s[i]=-s[i];}
        else {double p=exp(-dE*rkbt); if(p>=ran3(&seed)) {s[i]=-s[i];}}
        if(j>=teq) {int jj=j-teq+1; if((jj%mod)==0) mM+=fabs(M());}
    }
    mM*=(rmts/double(mod));
}

void ising2d::SMS(double t, double b)
{
    T=t; B=b;
    int seed=-845346743;
    double mss,uss,rkbt=1./(kB*T);
    bool flip=false;
    for(int j=0;j<ns;j++)
    {
        int i=int(double((N-1))*ran3(&seed));
        double dE=2.*double(s[i])*((J*SumNN(i))+(B*mo));
        if(dE<=0) {s[i]=-s[i]; flip=true;}
        else {double p=exp(-dE*rkbt); if(p>=ran3(&seed)) {s[i]=-s[i]; flip=true;}}
        if(j==teq) {initM(mss); initU(uss);}
        if(j>teq) {if(flip) {update1M(mss,i); update1U(uss,dE); flip=false;} update2M(mss); update2U(uss);}
    }
    rescale();
}

double ising2d::SumNN(int i)
{
    double sum=0;
    for(int j=0;j<Nbrs;j++) sum+=double(s[nn[i][j]]);
    return sum;
}

double ising2d::SumSpinProduct() {double sum=0; for(int i=0;i<N;i++) {sum+=double(s[i])*SumNN(i);} return sum;}

double ising2d::SumSpins() {double sum=0; for(int i=0;i<N;i++) {sum+=double(s[i]);} return sum;}

double ising2d::U() {return (-(0.5*J*SumSpinProduct())-(B*mo*SumSpins()));}

double ising2d::u() {return (U()*rN);}

double ising2d::M() {return (SumSpins()*rN);}

void ising2d::rescale() {mM*=rmts; mU*=rmts;}

double ising2d::avgM() {return mM;}

double ising2d::avgU() {return mU;}

double ising2d::avgu() {return (avgU()*rN);}

double ising2d::avgC() {return (((mU2*rmts)-(pow(mU,2)))/(kB*pow(T,2)));}

double ising2d::avgc() {return (avgC()*rN);}

void ising2d::initM(double &a) {mM=a=M(); /*mM=fabs(a);*/ mM2=pow(a,2);}

void ising2d::update1M(double &a, int i) {a+=(2.*double(s[i])*rN);}

void ising2d::update2M(double a) {mM+=a;/*mM+=fabs(a);*/ mM2+=pow(a,2);}

void ising2d::initU(double &a) {mU=a=U(); mU2=pow(a,2);}

void ising2d::update1U(double &a, double e) {a+=e;}

void ising2d::update2U(double a) {mU+=a; mU2+=pow(a,2);}

double ising2d::avgX() {return (((mM2*rmts)-(pow(mM,2)))/(kB*T));}



double ising2d::ran3(int *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;
    
	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}