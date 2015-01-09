/** @file evap.c
	@brief Evaporation routines.
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>


#include "mol.h"
/**Evaporate molecules from the system,  this functions takes molecules from the edge of the box
@arg in The system to be acted upon
@arg name the name of the solvent molecule
@arg the proportion of the molecule to be removed 1.0-0.0
*/
void evap_edge(struct system * in,char *name,double percent)
{
if (percent==0) return;
//first create list of min distances
int i;
int pos=0;
double *x=NULL;
double *y=NULL;
double *z=NULL;

int len=0;

for (i=0;i<in->length;i++)
{
	if (strcmp(in->mols[i].resid,name)==0)
	{
		len++;
	}
}

if (len==0)
{
printf("No solvent found in file\n");
return;
}

int remove=len*percent;
if ((remove==0)&&(len==1)) remove=1;

x=(double*)realloc(x,in->length*sizeof(double));
y=(double*)realloc(y,in->length*sizeof(double));
z=(double*)realloc(z,in->length*sizeof(double));
size_t * px = malloc (in->length * sizeof(size_t));
size_t * py = malloc (in->length * sizeof(size_t));
size_t * pz = malloc (in->length * sizeof(size_t));
pos=0;

for (i=0;i<in->length;i++)
{
		x[pos]=in->mols[i].cog.x;
		y[pos]=in->mols[i].cog.y;
		z[pos]=in->mols[i].cog.z;
		pos++;
}


gsl_sort_index (px, x, 1, (size_t)in->length);
gsl_sort_index (py, y, 1, (size_t)in->length);
gsl_sort_index (pz, z, 1, (size_t)in->length);




int removed=0;
int xpos=0;
int ypos=0;
int zpos=0;
int mol=0;
int cutoff=0;
	do
	{

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
		cutoff=0;
		if ((xpos<in->length)&&((in->length-1-xpos)>=0))
		{
			mol=px[xpos];
			if (strcmp(in->mols[mol].resid,name)==0)
			{
				if (in->mols[mol].del==FALSE)
				{
					in->mols[mol].del=TRUE;
					removed++;
				}
			}

			if (removed>=remove) break;

			mol=px[in->length-1-xpos];
			if (strcmp(in->mols[mol].resid,name)==0)
			{
				if (in->mols[mol].del==FALSE)
				{
					in->mols[mol].del=TRUE;
					removed++;
				}
			}

			if (removed>=remove) break;
			xpos++;
		}else
		{
			cutoff++;
		}		
		

//yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy

		if ((ypos<in->length)&&((in->length-1-ypos)>=0))
		{
			mol=py[ypos];
			if (strcmp(in->mols[mol].resid,name)==0)
			{
				if (in->mols[mol].del==FALSE)
				{
					in->mols[mol].del=TRUE;
					removed++;
				}
			}

			if (removed>=remove) break;

			mol=py[in->length-1-ypos];
			if (strcmp(in->mols[mol].resid,name)==0)
			{
				if (in->mols[mol].del==FALSE)
				{
					in->mols[mol].del=TRUE;
					removed++;
				}
			}

			if (removed>=remove) break;
			ypos++;
		}else
		{
			cutoff++;
		}

//zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

		if ((zpos<in->length)&&((in->length-1-zpos)>=0))
		{
			mol=pz[zpos];
			if (strcmp(in->mols[mol].resid,name)==0)
			{
				if (in->mols[mol].del==FALSE)
				{
					in->mols[mol].del=TRUE;
					removed++;
				}
			}

			if (removed>=remove) break;

			mol=pz[in->length-1-zpos];
			if (strcmp(in->mols[mol].resid,name)==0)
			{
				if (in->mols[mol].del==FALSE)
				{
					in->mols[mol].del=TRUE;
					removed++;
				}
			}

			if (removed>=remove) break;
			zpos++;	
		}else
		{
			cutoff++;
		}

	}while((removed<remove)&&(cutoff<3));


free(x);
free(y);
free(z);

free(px);
free(py);
free(pz);
}


/** Evaporate molecules which are the furthest away from the other molecules
@arg in The system to be acted upon
@arg name the name of the solvent molecule
@arg the proportion of the molecule to be removed 1.0-0.0
*/
void evap_min(struct system * in,char *name,double percent)
{
if (percent==0) return;
//first create list of min distances
int i;
int ii;
double *min=NULL;
int len=0;
struct vec delta;
double dist=10000.0;
printf("here\n");
for (i=0;i<in->length;i++)
{
	if (strcmp(in->mols[i].resid,name)==0)
	{
		len++;
	}
}
printf("here\n");
min=(double*)realloc(min,in->length*sizeof(double));
for (i=0;i<in->length;i++)
{
min[i]=10000;
}
printf("here\n");
for (i=0;i<in->length;i++)
{

	for (ii=0;ii<in->length;ii++)
	{
		cpy_vec(&delta,&(in->mols[i].cog));
		sub_vec(&delta,&(in->mols[ii].cog));
		dist=mod_vec(&delta);
		if (ii!=i)
		if (dist<min[i])
		{
			min[i]=dist;
			//printf("%lf %lf\n",dist,min[i]);
			//getchar();
		}
	}
//	min[i]=get_random_number(100);
	//printf("%lf %lf\n",dist,min[i]);
	//printf("%lf %lf\n",dist,min[i]);
	if ((i%100)==0) printf("%d/%d %lf\n",i,in->length,((double)i/(double)in->length)*100.0);

}
printf("here\n");

size_t * p = malloc (in->length * sizeof(size_t));

gsl_sort_index (p, min, 1, (size_t)in->length);
printf("herea\n");
if (len==0)
{
printf("No solvent found in file\n");
}
else
{
	int remove=len*percent;
	if ((remove==0)&&(len==1)) remove=1;
	int removed=0;
	int pos=in->length-1;
	int mol=0;

	do
	{
		mol=p[pos];
		if (strcmp(in->mols[mol].resid,name)==0)
		{
			if (in->mols[mol].del==FALSE)
			{
				in->mols[mol].del=TRUE;
				removed++;
			}
		}
	pos--;

	}while((removed<remove)&&(pos>0));
}

}

/**Evaporate molecules from the system,  molecules are removed randomly
@arg in The system to be acted upon
@arg name the name of the solvent molecule
@arg the proportion of the molecule to be removed 1.0-0.0
*/
void evap(struct system * in,char *name,double percent)
{
if (percent==0) return;
int i;
int *array=NULL;
int len=0;
for (i=0;i<in->length;i++)
{
	if (strcmp(in->mols[i].resid,name)==0)
	{
		len++;
		array=(int*)realloc(array,len*sizeof(int));
		array[len-1]=i;
	}
}

if (len==0)
{
printf("No solvent found in file\n");
}
else
{
	int remove=len*percent;
	if ((remove==0)&&(len==1)) remove=1;
	int removed=0;
	int rnd;
	do
	{
		rnd=(int)get_random_number(len);
		if (rnd==len) rnd=len-1;



		if (in->mols[array[rnd]].del==FALSE)
		{
			in->mols[array[rnd]].del=TRUE;
			removed++;
		}

	}while(removed<remove);
}

free(array);
}
