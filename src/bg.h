#ifndef BEDGRAPH_H
#define BEDGRAPH_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

struct bgdata
{
    std::vector<std::string> chr;
	std::vector<std::vector<int> > start;
    std::vector<std::vector<int> > end;
	std::vector<std::vector<float> > val;
    std::vector<float> v;

    int ch;
    int res;
    int len;
    int load(char *fn);
	bool set(std::string &c, int r);
	bool set(std::vector<float> &o, std::string &c, int r);
    void getlevel(std::vector<float> &o, std::string &c, int s, int e, int r, bool isaverage);
    void getcoverage(std::vector<float> &o, std::string &c, int s, int e, int r);
    int getChrID(std::string &c);
    float &operator[](int i);
    float operator()(int c,int i);
    float operator()(std::string &c,int i);
};

int bgdata::load(char *fn)
{
	std::ifstream in(fn);
    int n=0; 
    while(true)
    {
        std::string c;
        in>>c;
	    if(in.eof()) break;
    	while(c[0]=='#'||c=="browser"||c=="track")
        {
            std::getline(in,c);
            in>>c;
	        if(in.eof()) break;
        }
        ch=getChrID(c);
        if(ch==-1)
        {
            ch=chr.size();
            chr.push_back(c);
            start.resize(ch+1);
            end.resize(ch+1);
            val.resize(ch+1);
        }
        int s,e;
        float vl;
        in>>s>>e>>vl;
        if(vl!=0)
        {
            start[ch].push_back(s);
            end[ch].push_back(e);
            val[ch].push_back(vl);
        }
        n++;
    }
    in.close();
    return n;
};


bool bgdata::set(std::string &c, int r=1)
{
    int cc=getChrID(c);
    if(cc<0&&cc>=chr.size()) return false;
    res=r;
    ch=cc;
    len=end[ch].back();
	v.resize(len/res+1);
    v.assign(v.size(),0);
    for(int i=0;i<val[ch].size();i++)
        for(int pos=start[ch][i];pos<end[ch][i];pos++)
            v[pos/res]+=val[ch][i]/res;
    return true;
};

bool bgdata::set(std::vector<float> &o, std::string &c, int r=1)
{
    int cc=getChrID(c);
    if(cc<0&&cc>=chr.size()) return false;
    len=end[cc].back();
	o.resize(len/r+1);
    o.assign(v.size(),0);
    for(int i=0;i<val[cc].size();i++)
        for(int pos=start[cc][i];pos<end[cc][i];pos++)
            v[pos/r]+=val[cc][i]/r;
    return true;
};

int bgdata::getChrID(std::string &c)
{
    std::string d="chr";
    if(c.find("chr")==0) d=c.substr(3);
    else d+=c;
    for(int i=0;i<chr.size();i++)
    {
        if(c==chr[i]||d==chr[i]) return i;
    }
    return -1;
};


float &bgdata::operator[] (int i)
{
   int r=i/res;
   if(r<0) r=0;
   else if(r>=v.size()) r=v.size()-1;
   return v[r];
};

float bgdata::operator() (int c, int i)
{
    int pos=std::upper_bound(start[c].begin(),start[c].end(),i)-start[c].begin();
    if(pos>0) if(end[c][pos-1]>i) return val[c][pos-1];
    return 0;
};

float bgdata::operator() (std::string &c, int i)
{
    int cc=getChrID(c);
    if(cc<0) return 0;
    int pos=std::upper_bound(start[cc].begin(),start[cc].end(),i)-start[cc].begin();
    if(pos>0) if(end[cc][pos-1]>i) return val[cc][pos-1];
    return 0;
};

void bgdata::getlevel(std::vector<float> &o, std::string &c, int s, int e, int r=1, bool isaverage=false)
{
    int cc=getChrID(c);
    if(cc<0&&cc>=chr.size()) return;
    int size=(e-s)/r+1;
	o.clear();
	o.resize(size);
    int startpos=std::upper_bound(end[cc].begin(),end[cc].end(),s)-end[cc].begin();
    int endpos=std::upper_bound(start[cc].begin(),start[cc].end(),e)-start[cc].begin();
    if(startpos==start[cc].size()) return;
    if(start[cc][startpos]>=e) return;
    for(int i=start[cc][startpos];i<end[cc][startpos];i++) if(i>=s) o[(i-s)/r]+=val[cc][startpos];
    for(int i=startpos+1;i<endpos-1;i++)
        for(int j=start[cc][i];j<end[cc][i];j++) o[(j-s)/r]+=val[cc][i];
    for(int i=start[cc][endpos-1];i<end[cc][endpos-1];i++) if(i>=s&&i<e) o[(i-s)/r]+=val[cc][endpos-1];

    if(isaverage)
    {
        std::vector<int> count((e-s)/r+1);
        for(int i=e;i<s;i++) count[(i-s)/r]++;
        for(int i=0;i<size;i++) if(count[i]>0) o[i]/=count[i];
    }
}

void bgdata::getcoverage(std::vector<float> &o, std::string &c, int s, int e, int r=1)
{
    int cc=getChrID(c);
    if(cc<0) return;
    int size=(e-s)/r+1;
    o.clear();
	o.resize(size);
    int startpos=std::upper_bound(end[cc].begin(),end[cc].end(),s)-end[cc].begin();
    int endpos=std::upper_bound(start[cc].begin(),start[cc].end(),e)-start[cc].begin();
    if(startpos==start[cc].size()) return;
    if(start[cc][startpos]>=e) return;
    for(int i=start[cc][startpos];i<end[cc][startpos];i++) if(i>=s&&val[cc][startpos]!=0) o[(i-s)/r]++;
    for(int i=startpos+1;i<endpos-1;i++)
        for(int j=start[cc][i];j<end[cc][i];j++) if(val[cc][i]!=0) o[(j-s)/r]++;
    for(int i=start[cc][endpos-1];i<end[cc][endpos-1];i++) if(i>=s&&i<e&&val[cc][endpos-1]!=0) o[(i-s)/r]++;

    std::vector<int> count((e-s)/r+1);
    for(int i=e;i<s;i++) count[(i-s)/r]++;
    for(int i=0;i<size;i++) if(count[i]>0) o[i]/=count[i];
}
#endif
