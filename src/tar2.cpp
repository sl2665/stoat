#include "bg.h"
#include "hmm1.h"
#include <ctime>
#include <cstdlib>
#include <cstring>

using namespace std;

double drand()
{
    return (double)rand()/RAND_MAX+(double)rand()/RAND_MAX/RAND_MAX;
}

int buildc(int *v, int *p, int len, int bin, bgdata &r, bgdata &m, int pos, char strand)
{
	for(int i=0;i<len/bin;i++) { v[i]=0; p[i]=-1; }
	int posit=0;
    switch(strand)
    {
        case '+':
            for(int i=0;i<len;i++) if(pos+i<r.len) if(m.v[pos+i]<1) r.v[pos+i]=-1;
            for(int i=0;i<len;i++)
			{
				if(pos+i<r.len) if(r.v[pos+i]>=0)
				{
					if(r.v[pos+i]>0) v[posit/bin]=1;
					if(p[posit/bin]==-1) p[posit/bin]=pos+i;
					posit++;
				}
			}
            break;
        case '-':
            for(int i=0;i<len;i++) if(pos-i>=0) if(m.v[pos-i]<1) r.v[pos-i]=1;
            for(int i=0;i<len;i++)
			{
				if(pos-i>=0) if(r.v[pos-i]<=0)
				{
					if(r.v[pos-i]<0) v[posit/bin]=1;
					if(p[posit/bin]==-1) p[posit/bin]=pos-i;
					posit++;
				}
			}
            break;
    }
	return posit/bin+1;
}


int main(int argc, char *argv[])
{
    if(argc<3)
    {
        cout<<"arguments:\n-p : plus strand bedgraph\n-m : minus strand bedgraph"<<endl;
        cout<<"-mp : mappability bedgraph"<<endl;
        cout<<"-b bin size"<<endl;
        return 0;
    }

    int arit=1;
    char *pfn, *mfn, *ofn, *mpfn;
    char *cfn;

    int bin=1;
    while(arit<argc)
    {
        if(strcmp(argv[arit],"-o")==0) ofn=argv[arit+1];
        else if(strcmp(argv[arit],"-p")==0) pfn=argv[arit+1];
        else if(strcmp(argv[arit],"-m")==0) mfn=argv[arit+1];
        else if(strcmp(argv[arit],"-mp")==0) mpfn=argv[arit+1];
        else if(strcmp(argv[arit],"-b")==0) bin=atoi(argv[arit+1]);
        arit+=2;
    }

    bgdata pb,mb,mp;
    cerr<<pb.load(pfn)<<" lines of plus strand data loaded"<<endl;
    cerr<<mb.load(mfn)<<" lines of minus strand data loaded"<<endl;
    cerr<<mp.load(mpfn)<<" lines of mappability data loaded"<<endl;
	ofstream out(ofn);
    srand(time(NULL));

	int it=1;
	cerr<<pb.chr.size()<<" chromosomes detected"<<endl;
	cerr<<"HMM estimates"<<endl;
	cerr<<"chromosome\tEmission1\tEmission2\tTransition1(kb)\tTransition2(kb)"<<endl;
	
	for(int i = 0; i < pb.chr.size(); ++i)
    {
        if(!pb.set(pb.chr[i])) continue;
        if(!mb.set(pb.chr[i])) continue;
        if(!mp.set(pb.chr[i])) continue;

		int len=pb.len;
		int *v, *p, *path;
		v=new int[len/bin+1];
		p=new int[len/bin+1];
		path=new int[len/bin+1];
		int posit, start, nlen;
		// plus strand
		nlen=buildc(v,p,len,bin,pb,mp,0,'+');
		hmm1 h1(nlen);
		h1.init(v,0.999999-pow(0.999,bin),0.999999-pow(0.9,bin),0.00001*bin, 0.00001*bin);
		h1.iteration('f','f',1000);
		cerr<<pb.chr[i]<<" plus\t";
		cerr<<h1.p_emit(0)/bin*1000<<"\t"<<h1.p_emit(1)/bin*1000<<"\t";
		cerr<<bin/h1.p_transit(0,1)/1000<<"\t"<<bin/h1.p_transit(1,0)/1000<<endl;
		h1.viterbi(path);
		posit=0, start=-1;
		while(1)
		{
			while(path[posit]==0&&posit<nlen) posit++;
			if(posit>=nlen) break;
			start=posit;
			while(path[posit]==1&&posit<nlen) posit++;
			cout<<pb.chr[i]<<"\t"<<p[start]<<"\t"<<p[posit]<<"\ttar#"<<it<<"\t.\t+"<<endl;
			if(posit>=nlen) break;
			it++;
		}
		
		// minus strand
		nlen=buildc(v,p,len,bin,mb,mp,len-1,'-');
		hmm1 h2(nlen);
		h2.init(v,0.999999-pow(0.999,bin),0.999999-pow(0.9,bin),0.00001*bin, 0.00001*bin);
		h2.iteration('f','f',1000);
		cerr<<pb.chr[i]<<" minus\t";
		cerr<<h2.p_emit(0)/bin*1000<<"\t"<<h2.p_emit(1)/bin*1000<<"\t";
		cerr<<bin/h2.p_transit(0,1)/1000<<"\t"<<bin/h2.p_transit(1,0)/1000<<endl;
		h2.viterbi(path);
		posit=0, start=-1;
		while(1)
		{
			while(path[posit]==0&&posit<nlen) posit++;
			if(posit>=nlen) break;
			start=posit;
			while(path[posit]==1&&posit<nlen) posit++;
			cout<<pb.chr[i]<<"\t"<<p[posit]<<"\t"<<p[start]<<"\ttar#"<<it<<"\t.\t-"<<endl;
			if(posit>=nlen) break;
			it++;
		}
		delete[] v;
	    delete[] p;
        delete[] path;
    }
    return 1;
}

