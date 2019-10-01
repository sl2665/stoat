#ifndef SMOOTH_H
#define SMOOTH_H

#include <cmath>
#include <algorithm>
#define _N_GAUSSIAN_TABLE 10

namespace smooth
{
    double *_gaussian_distribution[_N_GAUSSIAN_TABLE];
    double _gaussian_std;
    int _gaussian_table_range;

    double polynomial_smooth_function(double x, double w)
	{
		double r=x/w;
		r*=r;
		r=1-r;
		r*=r*r*36/35/w;
		return r;
	};
	
    template <class classType>
    void set_smooth(classType w, int table_id=0)
    {
        _gaussian_std=(double)w/2;
        _gaussian_table_range=2*_gaussian_std;
        _gaussian_distribution[table_id]=new double[2*_gaussian_table_range+1];

        if(_gaussian_table_range==0)
        {
            _gaussian_distribution[table_id][0]=1;
            return;
        }

        double x,s=0;
        for(int i=-_gaussian_table_range;i<=_gaussian_table_range;i++)
        {
            x=(double)i;
            x=x*x/_gaussian_std/_gaussian_std/2;
            x=exp(-x);
            s+=x;
            _gaussian_distribution[table_id][i+_gaussian_table_range]=x;
        }

                                    
        for(int i=-_gaussian_table_range;i<=_gaussian_table_range;i++)
            _gaussian_distribution[table_id][i+_gaussian_table_range]/=s;
    };
	
    template <class classType>
    void set_polynomial_smooth(classType w, int table_id=0)
    {
        _gaussian_table_range=(int)w;
        _gaussian_distribution[table_id]=new double[2*_gaussian_table_range+1];
		
        if(_gaussian_table_range==0)
        {
            _gaussian_distribution[table_id][0]=1;
            return;
        }
		
        double x,s=0;
        for(int i=-_gaussian_table_range;i<=_gaussian_table_range;i++)
        {
            x=polynomial_smooth_function((double)i,(double)w);
            s+=x;
            _gaussian_distribution[table_id][i+_gaussian_table_range]=x;
        }
		
        for(int i=-_gaussian_table_range;i<=_gaussian_table_range;i++)
            _gaussian_distribution[table_id][i+_gaussian_table_range]/=s;
    };	

    template<class classType>
    void add_smooth(classType *array, classType value, int pos, int lim, int table_id=0)
    {
        double *dstr=_gaussian_distribution[table_id];
        for(int i=-_gaussian_table_range, p=pos-_gaussian_table_range;i<=_gaussian_table_range;i++,p++)
            if(p>=0&&p<lim) array[p]+=dstr[i+_gaussian_table_range]*value;
    };

    template<class classType1, class classType2>
    void smooth_to(classType1 *dest, classType2 *source, int lim, int table_id=0)
    {
        double *dstr=_gaussian_distribution[table_id];
        for(int i=0;i<lim;i++) dest[i]=0;
        for(int i=0;i<lim;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[i];
        }
        for(int i=-_gaussian_table_range;i<0;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[0];
        }
        for(int i=lim;i<=lim+_gaussian_table_range;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[lim-1];
        }
    };

    template<class classType>
    void smooth_it(classType *source, int lim, int table_id=0)
    {
        double *dstr=_gaussian_distribution[table_id];
        classType *dest;
        dest=new classType[lim];
        std::fill_n(dest,lim,0);
        for(int i=0;i<lim;i++)
        {
            if(source[i]!=0) for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[i];
        }
        if(source[0]!=0) for(int i=-_gaussian_table_range;i<0;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[0];
        }
        if(source[lim-1]!=0) for(int i=lim;i<=lim+_gaussian_table_range;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[lim-1];
        }
        copy(dest,dest+lim,source);
        delete[] dest;
    };

    template<class classType>
    void smooth_it(std::vector<classType> &source, int table_id=0)
    {
        int lim=source.size();
        double *dstr=_gaussian_distribution[table_id];
        classType *dest;
        dest=new classType[lim];
        std::fill_n(dest,lim,0);
        for(int i=0;i<lim;i++)
        {
            if(source[i]!=0) for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[i];
        }
        if(source[0]!=0) for(int i=-_gaussian_table_range;i<0;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[0];
        }
        if(source[lim-1]!=0) for(int i=lim;i<=lim+_gaussian_table_range;i++)
        {
            for(int j=-_gaussian_table_range, p=i-_gaussian_table_range;j<=_gaussian_table_range;j++,p++)
                if(p>=0&&p<lim) dest[p]+=dstr[j+_gaussian_table_range]*source[lim-1];
        }
        copy(dest,dest+lim,source.begin());
        delete[] dest;
    };
}
#endif
