#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include <iostream>
using namespace std;

#include "BCTool.h"
#include "CStemCell.h"
#include "System.h"
#include "Random.h"

struct IMD _MD;
CRandom Rand;

void ReadIPF(char *fn);
void SetParVal(char *str, char const *conststr, char val[]);
void help();
void SetDefault();
void OutputParameter();

int main(int argc, char *argv[])
{
	CSystem *sys;
  	if (argc<2)
  	{
  		help();
		exit(0);
  	}
	
	ReadIPF(argv[1]);
	if(_MD.seed > 0) Rand.Initialized(_MD.seed);

    sys = new CSystem(_MD.N0);
    if(sys->Initialized())
    {
        sys->Run();
        return(1);
    }
    else
    {
        return(0);
    }
}

void help()
{
	cout<<"Usage: "<<endl;
	cout<<"bct_StemCell inputfile"<<endl;
	cout<<"inputfile: The name of input file."<<endl;
	cout<<" For example: bct md.in."<<endl;
}

void ReadIPF(char *fn)
{
	FILE *fp;
	char str[StrLength], *pst;
	if((fp = fopen(fn,"r"))==NULL)
	{
		cout<<"Cannot open the input file."<<endl;
		exit(0);
	}
	SetDefault();
	rewind(fp);
	while(!feof(fp))
	{
		fgets(str,StrLength,fp);
		if(str[0]=='#'){ continue;}
		SetParVal(str, "mdcrd=\"", _MD.mdcrd);
		SetParVal(str, "cellpar=\"",_MD.cellpar);
        if((pst=strstr(str,"dt="))!=NULL)
		{
			_MD.dt=atof(pst+3);
		}
		if((pst=strstr(str,"T1="))!=NULL)
		{
			_MD.T1=atof(pst+3);
		}
        if((pst=strstr(str,"T0="))!=NULL)
        {
            _MD.T0=atof(pst+3);
        }
        if((pst=strstr(str,"ntpx="))!=NULL)
		{
            _MD.ntpx=atoi(pst+5);
        }
        if((pst=strstr(str,"ntpr="))!=NULL)
		{
            _MD.ntpr=atoi(pst+5);
        }
		if((pst=strstr(str,"seed="))!=NULL)
		{
			_MD.seed=atoi(pst+5);
		}
        if((pst=strstr(str,"N0="))!=NULL)
        {
            _MD.N0=atoi(pst+3);
        }
	}
	fclose(fp);
}

void SetParVal(char *str, char const *conststr, char val[])
{
	char *pst;
	if((pst=strstr(str,conststr))!=NULL)
	{
		strcpy(val,pst+strlen(conststr));
		if((pst = strstr(val,"\""))!=NULL)
		{
			val[pst - val] = '\0';
		}			
	}
	return;
}

void OutputParameter()
{
//     cout<<"mdcrd="<<_MD.mdcrd<<endl;
//	 cout<<"N="<<_MD.N<<endl;
//	 cout<<"dt="<<_MD.dt<<endl;
//     cout<<"T1="<<_MD.T1<<endl;
//	 cout<<"ntpr="<<_MD.ntpr<<endl;
//	 cout<<"ntpx="<<_MD.ntpx<<endl;
}

void SetDefault()
{
     _MD.ntpx=0;
     _MD.ntpr=0;
     _MD.N0=1;
	 _MD.seed=0;
}
