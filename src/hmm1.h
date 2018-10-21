#include <iostream>
#include <fstream>
#include <cmath>

#define NSTATE 2
#define INF 1.0e64
#define NOBS 2
using namespace std;

// 2 state hmm
// binary observation values
// Parameter estimation using Baum-Welch algorithm

class hmm1
{
    // initial state probabilities
    double initial[NSTATE];
    // transition probabilities
    double transition[NSTATE][NSTATE];
    // emission probabilities
    double emission[NSTATE][NOBS];
    // number of base positions
    int n;
    // relative density
    int *A;

    inline double start(int s) { return initial[s]; };
    inline double go(int t, int s) { return transition[t][s]; };
    inline double out(int s, int a) { return emission[s][a]; };

    double *alphamatrix[NSTATE];
    double *betamatrix[NSTATE];
    double *gammamatrix[NSTATE][NSTATE];
    double *deltamatrix[NSTATE];

    double PrA;

    inline double &alpha(int j, int s) { return alphamatrix[s][j]; };
    inline double &beta(int j, int s) { return betamatrix[s][j]; };
    inline double &gamma(int j, int s, int t) { return gammamatrix[s][t][j]; };
    inline double &delta(int j, int s) { return deltamatrix[s][j]; };

    void getalpha();
    void getbeta();
    void getgamma();
    void getdelta();
    void nextHMM();
	void boundary();
   
public:
	hmm1(int Nbases);
    void init(int *v, double em0, double em1, double tr01, double tr10);		// initialize reads ratio relative to reference
    double logL();
    int iteration(char verb, char bound,int round);
    void viterbi(int *path);
	double p_emit(int s) { return exp(emission[s][1]); };
	double p_transit(int t, int s) { return transition[t][s]; };
    ~hmm1();
};

hmm1::hmm1(int Nbases)
{
    n=Nbases;
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
}

void hmm1::init(int *v, double em0=0.001, double em1=0.1, double tr01=0.0001, double tr10=0.001)
{
    for(int i=0;i<n;i++)
	{
		A[i]=v[i];
		if(A[i]>1) A[i]=1;
	}
	initial[0]=0.999;
	initial[1]=0.001;
	emission[0][0]=log(1.0-em0);
	emission[0][1]=log(em0);
	emission[1][0]=log(1.0-em1);
	emission[1][1]=log(em1);
	transition[0][0]=1.0-tr01;
	transition[0][1]=tr01;
	transition[1][0]=tr10;
	transition[1][1]=1.0-tr10;
};

void hmm1::boundary()
{
};

hmm1::~hmm1()
{
    delete[] A;
    for(int i=0;i<NSTATE;i++)
    {
        delete[] alphamatrix[i];
        delete[] betamatrix[i];
        for(int j=0;j<NSTATE;j++) delete[] gammamatrix[i][j];
    }
};

void hmm1::getalpha()
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

void hmm1::getbeta()
{
    double b[NSTATE],bv,sum;
    for(int s=0;s<NSTATE;s++) beta(n-1,s)=0;

    for(int j=n-2;j>=0;j--)
    {
        for(int s=0;s<NSTATE;s++)
        {
            bv=-INF;
            for(int u=0;u<NSTATE;u++)
			{
                if(go(s,u)>0) b[u]=log(go(s,u))+out(u,A[j+1])+beta(j+1,u);
                else b[u]=-INF;
				if(b[u]>bv) bv=b[u];
            }
            sum=0;
            for(int u=0;u<NSTATE;u++) sum+=exp(b[u]-bv);
            beta(j,s)=log(sum)+bv;
        }
    }
};

void hmm1::getgamma()
{
    for(int j=0;j<n-1;j++)
        for(int s=0;s<NSTATE;s++)
            for(int t=0;t<NSTATE;t++)
                if(go(s,t)>0) gamma(j,s,t)=alpha(j,s)+log(go(s,t))+out(t,A[j+1])+beta(j+1,t)-PrA;
                else gamma(j,s,t)=-INF;
};

void hmm1::getdelta()
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

double hmm1::logL()
{
    return PrA/log(10);
};

void hmm1::nextHMM()
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

        for(int a=0;a<NOBS;a++)
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


int hmm1::iteration(char verb='f',char bon='f',int rdl=200)
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
        if(verb=='t'||verb=='T')
		{
			cout<<"Iteration round "<<rd<<"\t log likelyhood "<<ll<<endl;
		}
		nextHMM();
    }
    
    return rd;
};


void hmm1::viterbi(int *path)
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
