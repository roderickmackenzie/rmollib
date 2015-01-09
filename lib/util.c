/** @file util.c
	@brief Basic utility functions
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include "mol.h"

void platofrm_test()
{
printf("struct comp_header=%d\n",sizeof (struct comp_header ));
printf("unsigned long=%d\n",sizeof (unsigned long ));
printf("uint64_t=%d\n",sizeof(uint64_t));
printf("double=%d\n",sizeof (double));
printf("int=%d\n",sizeof (int));
printf("char=%d\n",sizeof (char));
printf("float=%d\n",sizeof (float));

}
/**
Returns the RAM of an attom
*/
double ret_ram(char *in)
{
	//printf("%s\n",in);
		if (strcmp(in,"H")==0)
		{
			return 1;
		}else
		if (strcmp(in,"C")==0)
		{
			return 12;
		}else
		if (strcmp(in,"CT1")==0)
		{
			return 12;
		}else
		if (strcmp(in,"Cl2")==0)
		{
			return 35.45300;	
		}else
		if (strcmp(in,"Cl3")==0)
		{
			return 35.45300;	
		}else
		if (strcmp(in,"Cl4")==0)
		{
			return 35.45300;	
		}else
		{
			//return 12.0;
			printf("Warning: Element %s not found - using carbon\n",in);
			exit(0);
		}
}



int strpartcmp(char *in, char *match)
{
int i;
if (strlen(in)<strlen(match)) return -1;
for (i=0;i<strlen(match);i++)
{
if (in[i]!=match[i]) return -1;
}
return 0;
}


double get_max(double *in,int length)
{
int i;
double max=in[0];
for (i=0;i<length;i++)
{
	if (in[i]>max)
	{
		max=in[i];
		//printf("%d %lf %lf\n",i,max,in[i]);
		//getchar();
	}
}
return max;
}

double round(double number)
{
char temp[100];
double out;

sprintf(temp,"%.3lf",number);
sscanf(temp,"%lf",&out);

return out;
}

int mycmp4(char *buf,char *in)
{
//printf("size = %d %d %d %s----%s\n",sizeof(int),*(int*)(buf),*((int*)in),buf,in);

if (strlen(buf)<4) return FALSE;
if ((*((int*)in)-*(int*)(buf))!=0) return FALSE;
return TRUE;
}

int mycmp(char *buf,char *in)
{
char *tin=in;
char *tbuf=buf;
while (!(*(unsigned char *) tin - *(unsigned char *) tbuf) && *tbuf)
{
tin++;
tbuf++;
}
if (*tin==0) return TRUE;

return FALSE;
}

void remove_white_space(char *in)
{
int i;
int ii=0;
for (i=0;i<strlen(in);i++)
{

if (in[i]!=' ')
{
	in[ii]=in[i];
	ii++;	
}

}
in[ii]=0;
}

///Print the copyright
void copyright()
{
printf("Copyright Roderick MacKenzie 2009 released under GPL v2. see www.gnu.org for license details\n ");

}

int scanarg( char* in[],int count,char * find)
{
int i;
for (i=0;i<count;i++)
{
if (strcmp(in[i],find)==0) return TRUE;
}
return FALSE;
}

int get_arg_plusone_pos( char* in[],int count,char * find)
{
int i;
for (i=0;i<count;i++)
{
if (strcmp(in[i],find)==0)
{
       if ((i+1)<count)
       {
               return i+1;
       }else
       {
               return -1;
       }
}
}
return -1;
}


char * get_arg_plusone( char* in[],int count,char * find)
{
int i;
static char no[] = "";
for (i=0;i<count;i++)
{
//printf("%s %s\n",in[i],find);
if (strcmp(in[i],find)==0)
{
       if ((i+1)<count)
       {
		//printf("%s\n",in[i+1]);
               return in[i+1];
       }else
       {
               return no;
       }
}
}

return no;
}

void right_align(char *out, int space,char *in)
{
	int len=strlen(in);
	int white=space-len;
	int i=0;
	for (i=0;i<white;i++)
	{
		strcat(out, " ");
	}
	strcat(out,in);
}



