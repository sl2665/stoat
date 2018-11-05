#include "bg.h"
#include "hmm5.h"
#include "bed.h"
#include <ctime>
#include <cstdlib>
#include <cstring>

using namespace std;

double drand()
{
    return (double)rand()/RAND_MAX+(double)rand()/RAND_MAX/RAND_MAX;
}

int buildc(float *v, int *p, int len, int bin, bgdata &r, bgdata &m, int pos, char strand)
{
	int *mp;
	mp = new int[len/bin+1];
	for(int i=0;i<len/bin+1;i++) { v[i]=0; p[i]=-1; mp[i]=0; }
	int posit=0;
	float average=0;
    switch(strand)
    {
        case '+':
            for(int i=0;i<len;i++) if(pos+i<r.len&&pos+i>=0) if(m.v[pos+i]<1)
				for(int k=pos+i-50;k<pos+i+50;k++) if(k>=0&&k<r.len) r.v[k]=-1;
            for(int i=0;i<len;i++)
			{
				if(pos+i<r.len&&pos+i>=0) if(r.v[pos+i]>=0)
				{
					if(r.v[pos+i]>0) v[posit/bin]+=r.v[pos+i];
					if(p[posit/bin]==-1) p[posit/bin]=pos+i;
					mp[posit/bin]++;
					posit++;
				}
			}
            break;
        case '-':
            for(int i=0;i<len;i++) if(pos-i>=0&&pos-i<r.len) if(m.v[pos-i]<1)
				for(int k=pos-i-50;k<pos-i+50;k++) if(k>=0&&k<r.len) r.v[k]=1;
            for(int i=0;i<len;i++)
			{
				if(pos-i>=0&&pos-i<r.len) if(r.v[pos-i]<=0)
				{
					if(r.v[pos-i]<0) v[posit/bin]-=r.v[pos-i];
					if(p[posit/bin]==-1) p[posit/bin]=pos-i;
					mp[posit/bin]++;
					posit++;
				}
			}
            break;
    }
	int blen=(posit-1)/bin+1;
	for(int i=0;i<blen;i++) if(mp[i]>0) v[i]/=mp[i];	// reads per mappable base
	for(int i=0;i<blen;i++) average+=v[i]/blen;
	if(average>0) for(int i=0;i<blen;i++) v[i]/=average;
	delete[] mp;
	return blen;
}


int main(int argc, char *argv[])
{
    if(argc<3)
    {
        cout<<"arguments"<<endl;
		cout<<"-p : plus strand bedgraph"<<endl;
		cout<<"-m : minus strand bedgraph"<<endl;
        cout<<"-mp : mappability bedgraph"<<endl;
        cout<<"-g : gene list bed"<<endl;
		cout<<"-bn : bin numbers(default=100)"<<endl;
        cout<<"-bs : min bin size(default=1000)"<<endl;
		return 0;
    }

    int arit=1;
    char *pfn, *mfn, *ofn, *mpfn, *gfn;
    char *cfn;

    int nbin=50,sbin=1000;
	float drop_off=0.1;
    while(arit<argc)
    {
        if(strcmp(argv[arit],"-p")==0) pfn=argv[arit+1];
        else if(strcmp(argv[arit],"-m")==0) mfn=argv[arit+1];
        else if(strcmp(argv[arit],"-mp")==0) mpfn=argv[arit+1];
        else if(strcmp(argv[arit],"-g")==0) gfn=argv[arit+1];
        else if(strcmp(argv[arit],"-bn")==0) nbin=atoi(argv[arit+1]);
        else if(strcmp(argv[arit],"-bs")==0) sbin=atoi(argv[arit+1]);
        arit+=2;
    }

    bgdata pb,mb,mp;
    cerr<<pb.load(pfn)<<" lines of plus strand data loaded"<<endl;
    cerr<<mb.load(mfn)<<" lines of minus strand data loaded"<<endl;
    cerr<<mp.load(mpfn)<<" lines of mappability data loaded"<<endl;
    
	srand(time(NULL));
	
	bedFile g;
	cerr<<g.load(gfn)<<" lines of gene annotation bed data loaded"<<endl;
   
    int it=1;

	for(int i = 0; i < pb.chr.size(); ++i)
    {
		string chr=pb.chr[i];
		cerr<<"Processing "<<chr<<"..."<<endl;
        if(!pb.set(chr)) continue;
        if(!mb.set(chr)) continue;
        if(!mp.set(chr)) continue;

		int chrlen=pb.len;
		
		float *v;
		int *p, *path;
		int nstart, nend, nlen, szbin, nsbin, posit;
		int flag;
		for(int i=0;i<g.n;i++)
		{
			if(g.t[i].chr==chr)
			{
				nstart=(3*g.t[i].start-g.t[i].end)/2;
				nend=(3*g.t[i].end-g.t[i].start)/2;
				nlen=nend-nstart;
				nlen=(nlen>0)?nlen:-nlen;
				
				if(nlen<nbin*sbin) nlen=nbin*sbin;
				szbin=nlen/nbin;

				v=new float[nbin+1];
				p=new int[nbin+1];
				path=new int[nbin+1];
				if(g.t[i].strand=='+') nsbin=buildc(v,p,nlen,szbin,pb,mp,nstart,'+');
				else nsbin=buildc(v,p,nlen,szbin,mb,mp,nend,'-');
				
				flag=0;
				for(int j=0;j<nsbin;j++) if(v[j]>0) flag=1;
				if(nsbin>10&&flag>0)
				{
					hmm5 h(nsbin);
					h.init(v);
					h.iteration();
					h.viterbi(path);
				
					posit=nsbin/4;
					while(path[posit]==0&&posit<nsbin-1) posit++;
					while(path[posit]==1&&posit<nsbin-1) posit++;
					posit--;
					if(g.t[i].strand=='+') g.t[i].end=p[posit];
					else g.t[i].start=p[posit];
					cout<<g.t[i];
				}
				delete[] v;
				delete[] p;
				delete[] path;	
			}
		}
	}
}
