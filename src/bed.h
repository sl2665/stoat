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
        if(!std::getline(in, input)) return(false);
        std::stringstream ss(input);
        ss >> chr >> start >> end;
        ss >> name >> score >> strand;
        ss >> cdsStart >> cdsEnd >> RGB;
        ss >> exonCount >> exonSizes >> exonStarts;
        return(true);
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
