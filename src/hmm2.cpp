#include "bg.h"
#include "bed.h"
#include "hmm2.h"
#include <ctime>
#include <cstdlib>
#include <cstring>

using namespace std;

double drand()
{
    return (double)rand()/RAND_MAX+(double)rand()/RAND_MAX/RAND_MAX;
}

void buildrc(float *v, int len, int bin, bgdata &r, bgdata &s, int pos, char strand)
{
    int npos;
    int *rv;
    rv=new int[len/bin+1];
    for(int i=0;i<len/bin;i++) rv[i]=v[i]=0;

    switch(strand)
    {
        case '+':
            // reference density per bin
            for(int i=0;i<len;i++) if(r.v[pos+i]>0) rv[i/bin]++;
            // local density normalization per bin
            for(int i=0;i<len;i++) if(s.v[pos+i]>0) v[i/bin]++;
            break;
        case '-':
            // reference density per bin
            for(int i=0;i<len;i++) if(r.v[pos-i]<0) rv[i/bin]++;
            // local density normalization per bin
            for(int i=0;i<len;i++) if(s.v[pos-i]<0) v[i/bin]++;
            break;
    }
	for(int i=0;i<len/bin;i++)
	{
		if(rv[i]>5) v[i]/=rv[i];
		else if(i>0) v[i]=v[i-1];
		else v[i]=0.1;
	}
    delete[] rv;
}


int main(int argc, char *argv[])
{
    if(argc<3)
    {
        cerr<<"arguments:\n\t-p : plus strand bedgraph\n\t-m : minus strand bedgraph"<<endl;
        cerr<<"\t-p0 : plus strand reference\n\t-m0 : minus strand reference"<<endl;
        cerr<<"\t-g : genelist bed"<<endl;
        cerr<<"\t-bs : bin size (default = 1000)\n\t-bc : bin count (default = 50)"<<endl;
        return 0;
    }

    int arit=1;
    char *pfn, *mfn, *ofn, *p0fn, *m0fn;
    char *gfn;

    int bin=1000, bc = 50;
    while(arit<argc)
    {
        if(strcmp(argv[arit],"-o")==0) ofn=argv[arit+1];
        else if(strcmp(argv[arit],"-p")==0) pfn=argv[arit+1];
        else if(strcmp(argv[arit],"-m")==0) mfn=argv[arit+1];
        else if(strcmp(argv[arit],"-p0")==0) p0fn=argv[arit+1];
        else if(strcmp(argv[arit],"-m0")==0) m0fn=argv[arit+1];
        else if(strcmp(argv[arit],"-g")==0) gfn=argv[arit+1];
        else if(strcmp(argv[arit],"-bs")==0) bin=atoi(argv[arit+1]);
        else if(strcmp(argv[arit],"-bc")==0) bc=atoi(argv[arit+1]);
        arit+=2;
    }

	// Read bedgraph data
    bgdata pb,mb,p0b,m0b;
    cerr<<pb.load(pfn)<<" lines of plus strand data loaded"<<endl;
    cerr<<mb.load(mfn)<<" lines of minus strand data loaded"<<endl;
    cerr<<p0b.load(p0fn)<<" lines of plus strand reference data loaded"<<endl;
    cerr<<m0b.load(m0fn)<<" lines of minus strand reference data loaded"<<endl;
    
    srand(time(NULL));
	
	cerr<<pb.chr.size()<<" chromosomes detected"<<endl;
	cout<<"id\tlen\trounds\ttransition\tdensity1\tdensity2"<<endl;

	// Read gene annotation bed data
	bedFile g;
	cerr<<g.load(gfn)<<" lines of gene annotation bed data loaded"<<endl;

	for(int i = 0; i < pb.chr.size(); ++i)
    {
		string chr=pb.chr[i];
		cerr<<"Processing "<<chr<<"..."<<endl;
        if(!pb.set(chr)) continue;
        if(!mb.set(chr)) continue;
        if(!p0b.set(chr)) continue;
        if(!m0b.set(chr)) continue;
		
        for(int j = 0;j < g.n; ++j)
        {
            if(g.t[j].chr==chr)
            {
                if(g.t[j].strand=='+')
                {	
                    int len = g.t[j].end - g.t[j].start - bin;
                    if(len>bc*bin) len=bc*bin;
					if(len>10*bin)
					{
						hmm2 h(len/bin);
						float *v;
						v=new float[len/bin+1];
						buildrc(v,len,bin,p0b,pb,g.t[j].start+bin,'+');
						h.init(v);
						cout<<g.t[j].name<<"\t"<<len+bin<<"\t"<<h.iteration()<<"\t";
						cout<<bin*h.tpos(1)+bin<<"\t";
						cout<<h.density(0)<<"\t"<<h.density(1)<<endl;
						delete[] v;
					}
                }
                else
                {
                    int len=g.t[j].end-g.t[j].start-bin;
                    if(len>bc*bin) len=bc*bin;
					if(len>10*bin)
					{
						hmm2 h(len/bin);
						float *v;
						v=new float[len/bin+1];
						buildrc(v,len,bin,m0b,mb,g.t[j].end-bin,'-');
						h.init(v);
						cout<<g.t[j].name<<"\t"<<len+bin<<"\t"<<h.iteration()<<"\t";
						cout<<bin*h.tpos(1)+bin<<"\t";
						cout<<h.density(0)<<"\t"<<h.density(1)<<endl;
						delete[] v;
					}
                }
            }
        }
    }
    return 1;
}

