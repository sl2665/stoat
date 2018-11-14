#ifndef BED_H
#define BED_H

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

struct bedTrack
{
    std::string chr;
    int start, end;
    std::string name;
    std::string score;
    char strand;
    int cdsStart, cdsEnd;
    std::string RGB;
    int exonCount;
    std::string exonSizes;
    std::string exonStarts;
    
    int read(std::ifstream& in)
    {
        std::string input;
        if(in.eof()) return(false);
		if(!std::getline(in, input)) return(false);
        std::stringstream ss(input);
        ss >> chr >> start >> end;
		if(chr.empty()) return(false);
        ss >> name >> score >> strand;
        ss >> cdsStart >> cdsEnd >> RGB;
        ss >> exonCount >> exonSizes >> exonStarts;
        return(true);
    }

	friend std::ostream& operator<<(std::ostream& os, const bedTrack t)
	{
		os<<t.chr<<"\t"<<t.start<<"\t"<<t.end<<"\t";
		os<<t.name<<"\t"<<t.score<<"\t"<<t.strand<<std::endl;
		return os;
	}
};

struct bedFile
{
    std::vector <bedTrack> t;
	int n;
    int load(char *fn)
    {
        bedTrack r;
        n = 0;
        std::ifstream in(fn);
        while(r.read(in))
        {
            t.push_back(r);
            ++n;
        };
        return(n);
    }
};
#endif
