#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include "bed.h"
#include "bg.h"
using namespace std;

int main(int argc, char *argv[])
{
	if(argc<2)
	{
		cout<<"usage"<<endl;
		cout<<"-i\t: input bed"<<endl;
		cout<<"-p\t: plus strand bedgraph"<<endl;
		cout<<"-m\t: minus strand bedgraph"<<endl;
		cout<<"-win\t: window (default=1000)"<<endl;
		cout<<"-res\t: resolution (default=1)"<<endl;
		return 0;
	}

	char *ifn, *ofn, *cfn, *pfn, *mfn;
	int arit=1;
	int win=1000;
	int res=1;
	while(arit<argc)
	{
		if(strcmp(argv[arit],"-i")==0) { ifn=argv[arit+1]; arit+=2; }
		else if(strcmp(argv[arit],"-p")==0) { pfn=argv[arit+1]; arit+=2; }
		else if(strcmp(argv[arit],"-m")==0) { mfn=argv[arit+1]; arit+=2; }
		else if(strcmp(argv[arit],"-win")==0) { win=atoi(argv[arit+1]); arit+=2; }
		else if(strcmp(argv[arit],"-res")==0) { res=atoi(argv[arit+1]); arit+=2; }
		else arit++;
	}
	
	bedFile in;
	bgdata pl, mn;
	cerr<<in.load(ifn)<<" lines of input bed loaded"<<endl;
	cerr<<pl.load(pfn)<<" lines of plus strand bedgraph loaded"<<endl;
	cerr<<mn.load(mfn)<<" lines of minus strand bedgraph loaded"<<endl;

	for(int cit = 0; cit < pl.chr.size(); ++cit)
	{
		string chr=pl.chr[cit];
		cerr<<"Processing "<<chr<<"..."<<endl;
        if(!pl.set(chr, res)) continue;
        if(!mn.set(chr, res)) continue;
		for(int i=0;i<in.n;i++)
		{
			if(in.t[i].chr==chr)
			{
				float max_v=0;
				if(in.t[i].strand=='+')
				{
					int len=pl.v.size();
					stringstream buf;
					int max_p=in.t[i].start;
					for(int j=in.t[i].start;j<in.t[i].end;j+=res)
					{
						float tv=0;
						for(int k=(j-win/2)/res;k<(j+win/2)/res;k++) if(k>=0&&k<len) tv+=pl.v[k];	
						if(tv>max_v) { max_v=tv; max_p=j; }
					}
					buf<<max_p<<":"<<max_v;
					in.t[i].score = buf.str().c_str();
				}
				else
				{
					int len=mn.v.size();
					stringstream buf;
					int max_p=in.t[i].start;
					for(int j=in.t[i].start;j<in.t[i].end;j+=res)
					{
						float tv=0;
						for(int k=(j-win/2)/res;k<(j+win/2)/res;k++) if(k>=0&&k<len) tv-=mn.v[k];	
						if(tv>max_v) { max_v=tv; max_p=j; }
					}
					buf<<max_p<<":"<<max_v;
					in.t[i].score = buf.str().c_str();
				}
				cout<<in.t[i];
			}
		}
	}
}

