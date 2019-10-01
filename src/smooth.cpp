#include "bedgraph.h"			// Custom reading bedgraph files
#include "genelist.h"			// Genelist (bed) input
#include "arguments.h"			// Parsing command line arguments
#include "smooth.h"

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <functional>
using namespace std;

// Declarations
struct smParam						// Arguments
{
	char *fn;			// Bedgraph filenames
	int bin;			// Display bins
	int res;			// Resolution
	double mFactor;		// Multiply factor
	// Initialization
	smParam() { res=20; bin=1; mFactor=1;};
	bool get(int argc, char *argv[])	// Function to parse arguments
	{
		arguments a;
		a.push("-i","<input bedgraph>",fn);
		a.push("-res","Smoothing resolution (default = 20)",res,true);
		a.push("-bin","Display bins (default = 1)",bin,true);
		a.push("-mf","Multiply factor (default = 1)",mFactor,true);
		return(a.get(argc,argv));
	};
};

int main(int argc,char *argv[])
{
	// Parsing arguments
	smParam a;
	if(!a.get(argc,argv)) return 0;
	smooth::set_smooth((double)a.res/a.bin);
	
	// Load bedgraph data
	bgdata bg;
	bg.load(a.fn);	
	double totalReads = abs((double)bg.tc) / 1000000;	
	for(int i = 0; i < bg.chr.size(); ++i) {
		bg.set(bg.chr[i], a.bin);
		smooth::smooth_it(bg.v);
		for(int j = 0; j < bg.v.size()-1; ++j)
			if(bg.v[j]!=0) cout<<bg.chr[i]<<"\t"<<j*a.bin<<"\t"<<(j+1)*a.bin<<"\t"<<bg.v[j]/totalReads*a.mFactor*1000/a.bin<<endl;
	}
}
