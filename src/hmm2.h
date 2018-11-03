#include <iostream>
#include <fstream>
#include <cmath>

#define NSTATE 2
#define INF 1.0e64

using namespace std;

// 2 state hmm
// Parameter estimation using Baum-Welch algorithm

class hmm2
{
    // initial state probabilities
    double initial[NSTATE];
    // transition probabilities
    double transition[NSTATE][NSTATE];
    // emission probabilities
    double *emission[NSTATE];
    // number of base positions
    int n;
    // relative density
    int *A;
	// reference density
    int Dref;
	// maximum density
	int Dmax;

    inline double start(int s) { return initial[s]; };
    inline double go(int t, int s) { return transition[t][s]; };
    inline double out(int s, int a) { return emission[s][a]; };

    double *alphamatrix[NSTATE];
    double *betamatrix[NSTATE];
    double *gammamatrix[NSTATE][NSTATE];
    double *deltamatrix[NSTATE];

    double *nCk;
    double PrA;

    public:
    inline double &alpha(int j, int s) { return alphamatrix[s][j]; };
    inline double &beta(int j, int s) { return betamatrix[s][j]; };
    inline double &gamma(int j, int s, int t) { return gammamatrix[s][t][j]; };
    inline double &delta(int j, int s) { return deltamatrix[s][j]; };

    hmm2(int Nbases,int refdensity, int maxdensity);
    void init(float *v);		// initialize reads ratio relative to reference
    void getalpha();
    void getbeta();
    void getgamma();
    void getdelta();
    void nextHMM();
    double density(int s);
    double tpos(int s);
    double logL();
    int iteration(char bound,int round);
    void boundary();
    void viterbi(int *path);
    ~hmm2();
};

hmm2::hmm2(int Nbases,int refdensity=20, int maxdensity=40)
{
    n=Nbases;
    Dref=refdensity;
	Dmax=maxdensity;
    A=new int[n];
    for(int i=0;i<NSTATE;i++)
    {
        alphamatrix[i]=new double[n];
        betamatrix[i]=new double[n];
        for(int j=0;j<NSTATE;j++) gammamatrix[i][j]=new double[n];
        deltamatrix[i]=new double[n];
    }
    
    for(int i=0;i<NSTATE;i++)
        for(int j=0;j<NSTATE;j++)
            transition[i][j]=0;
    nCk=new double[Dmax+1];
    nCk[0]=0;
    for(int i=1;i<=Dmax;i++) nCk[i]=nCk[i-1]+log((float)(Dmax-i+1)/i);
	for(int i=0;i<NSTATE;i++) emission[i]=new double[Dmax+1];
}

void hmm2::init(float *v)
{
    for(int i=0;i<n;i++)
	{
		A[i]=v[i]*Dref;
		if(A[i]>Dmax) A[i]=Dmax;
	}
    initial[0]=1;
    initial[1]=0;
    transition[0][1]=0.02;   // transition 
    transition[0][0]=1.0-transition[0][1];
    transition[1][0]=0;
    transition[1][1]=1;

    double p[2];
	p[0]=0.1*Dref/Dmax;
	p[1]=1.0*Dref/Dmax;
    for(int s=0;s<NSTATE;s++) for(int a=0;a<=Dmax;a++) emission[s][a]=nCk[a]+log(p[s])*a+log(1-p[s])*(Dmax-a);
};

void hmm2::boundary()
{
    initial[0]=1;
    initial[1]=0;
    transition[1][1]=1;
    transition[1][0]=0;
    double p;
    for(int s=0;s<NSTATE;s++)
    {
        p=density(s)*Dref/Dmax;
        for(int a=0;a<=Dmax;a++) emission[s][a]=nCk[a]+log(p)*a+log(1-p)*(Dmax-a);
    }
};

hmm2::~hmm2()
{
    delete[] A;
    for(int i=0;i<NSTATE;i++)
    {
        delete[] alphamatrix[i];
        delete[] betamatrix[i];
        for(int j=0;j<NSTATE;j++) delete[] gammamatrix[i][j];
    	delete[] emission[i];
	}
};

void hmm2::getalpha()
{
    double a[NSTATE],av,sum;
    for(int s=0;s<NSTATE;s++) if(start(s)>0) alpha(0,s)=log(start(s))+out(s,A[0]); else alpha(0,s)=-INF;

    for(int j=1;j<n;j++)
    {
        for(int s=0;s<NSTATE;s++)
        {
            av=-INF;
			for(int t=0;t<NSTATE;t++)
            {
                if(go(t,s)>0) a[t]=alpha(j-1,t)+log(go(t,s))+out(s,A[j]);
                else a[t]=-INF;
				if(a[t]>av) av=a[t];
            }
            sum=0;
			for(int t=0;t<NSTATE;t++) sum+=exp(a[t]-av);
			alpha(j,s)=log(sum)+av;
        }
    }

    av=-INF;
    for(int s=0;s<NSTATE;s++) if(alpha(n-1,s)>av) av=alpha(n-1,s);
    sum=0;
    for(int s=0;s<NSTATE;s++) sum+=exp(alpha(n-1,s)-av);
    PrA=log(sum)+av;
};

void hmm2::getbeta()
{
    double b[NSTATE],bv,sum;
    for(int s=0;s<NSTATE;s++) beta(n-1,s)=0;

    for(int j=n-2;j>=0;j--)
    {
        for(int s=0;s<NSTATE;s++)
        {
            bv=-INF;
            for(int u=0;u<NSTATE;u++)
                if(go(s,u)>0)
                {
                    b[u]=log(go(s,u))+out(u,A[j+1])+beta(j+1,u);
                    if(b[u]>bv) bv=b[u];
                }
            sum=0;
            for(int u=0;u<NSTATE;u++) if(go(s,u)>0) sum+=exp(b[u]-bv);
            beta(j,s)=log(sum)+bv;
        }
    }
};

void hmm2::getgamma()
{
    for(int j=0;j<n-1;j++)
        for(int s=0;s<NSTATE;s++)
            for(int t=0;t<NSTATE;t++)
                if(go(s,t)>0) gamma(j,s,t)=alpha(j,s)+log(go(s,t))+out(t,A[j+1])+beta(j+1,t)-PrA;
                else gamma(j,s,t)=-INF;
};

void hmm2::getdelta()
{
    double d[NSTATE],dv,sum;
    for(int j=0;j<n-1;j++)
    {
        for(int s=0;s<NSTATE;s++)
        {
            dv=-INF;
            for(int u=0;u<NSTATE;u++)
            {
                d[u]=gamma(j,s,u);
                if(d[u]>dv) dv=d[u];
            }
            sum=0;
            for(int u=0;u<NSTATE;u++) sum+=exp(d[u]-dv);
            delta(j,s)=log(sum)+dv;
        }
    }
    for(int s=0;s<NSTATE;s++) delta(n-1,s)=alpha(n-1,s)-PrA;
};

double hmm2::logL()
{
    return PrA/log(10);
};

void hmm2::nextHMM()
{
    double Kv,sK,K,sum;
    Kv=-INF;
    for(int s=0;s<NSTATE;s++) if(delta(0,s)>Kv) Kv=delta(0,s);
    sum=0;
    for(int s=0;s<NSTATE;s++) sum+=exp(delta(0,s)-Kv);
    K=log(sum)+Kv;
    for(int s=0;s<NSTATE;s++) initial[s]=exp(delta(0,s)-K);   
    
    for(int s=0;s<NSTATE;s++)
    {
        Kv=-INF;
        for(int j=0;j<n-1;j++)
            for(int t=0;t<NSTATE;t++)
                if(gamma(j,s,t)>Kv) Kv=gamma(j,s,t);
        sum=0;
        for(int j=0;j<n-1;j++)
            for(int t=0;t<NSTATE;t++) sum+=exp(gamma(j,s,t)-Kv);
        sK=log(sum)+Kv;

        for(int t=0;t<NSTATE;t++)
        {
            Kv=-INF;
            for(int j=0;j<n-1;j++)
                if(gamma(j,s,t)>Kv) Kv=gamma(j,s,t);
            sum=0;
            for(int j=0;j<n-1;j++) sum+=exp(gamma(j,s,t)-Kv);
            K=log(sum)+Kv;
            
            transition[s][t]=exp(K-sK);
        }
    }

    for(int s=0;s<NSTATE;s++)
    {
        Kv=-INF;
        for(int j=0;j<n;j++)
            if(delta(j,s)>Kv) Kv=delta(j,s);
        sum=0;
        for(int j=0;j<n;j++) sum+=exp(delta(j,s)-Kv);
        sK=log(sum)+Kv;

        for(int a=0;a<=Dmax;a++)
        {
            Kv=-INF;
            for(int j=0;j<n;j++) if(A[j]==a) if(delta(j,s)>Kv) Kv=delta(j,s);
            sum=0;
            for(int j=0;j<n;j++) if(A[j]==a) sum+=exp(delta(j,s)-Kv);
            if(sum>0) K=log(sum)+Kv;
			else K=-INF;
            emission[s][a]=K-sK;
        } 
    }
};


int hmm2::iteration(char bon='f',int rdl=200)
{
    double ll=-0.9*INF, llpre=-INF;
    int it;
    int rd=0;

    while(ll>llpre&&rd<rdl)
    {
        rd++;
        if(bon=='t'||bon=='T') boundary();
        getalpha();
        getbeta();
        getgamma();
        getdelta();
        llpre=ll;
        ll=logL();
		nextHMM();
    }
    
    return rd;
};

double hmm2::density(int s)
{
    double d=0;
    for(int a=0;a<=Dmax;a++) d+=exp(out(s,a))*a/Dref;
    return d;
};

double hmm2::tpos(int s)
{
	return 1.0/go(s-1,s);
};

void hmm2::viterbi(int *path)
{
	double *vp[NSTATE];
	
	for(int s=0;s<NSTATE;s++)
	{
		vp[s]=new double[n];
		if(start(s)>0) vp[s][0]=out(s,A[0])+log(start(s));
		else vp[s][0]=-INF;
	}
	for(int j=1;j<n;j++)
	{
		for(int s=0;s<NSTATE;s++)
		{
			double maxvp=-INF;
			for(int u=0;u<NSTATE;u++)
			{
				if(go(u,s)>0) if(vp[u][j-1]+log(go(u,s))>maxvp) maxvp=vp[u][j-1]+log(go(u,s));
			}
			vp[s][j]=maxvp+out(s,A[j]);
		}
	}
	for(int j=0;j<n;j++)
	{
		path[j]=0;
		for(int s=1;s<NSTATE;s++) if(vp[s][j]>vp[path[j]][j]) path[j]=s;
	}
	for(int s=0;s<NSTATE;s++) delete[] vp[s];
};
