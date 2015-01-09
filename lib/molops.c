/** @file molops.c
	@brief Functions to operate directly on molecules
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include "mol.h"

/**Check the bonds assosiated with a molecule
@param my_mol1 molecule 1
@param my_mol2 molecule 2
@param chum1 dont know
@param chum2 dont know
@return dont know
*/
int mol_check_bond(struct mol *my_mol1,int chum1,struct mol *my_mol2,int chum2)
{
int i;
int found=FALSE;
for (i=0;i<my_mol1->nchums;i++)
{
	if (my_mol1->chums[i]==chum2)
	{
		found=TRUE;
		return TRUE;
	}
}

my_mol1->nchums++;
my_mol1->chums= (int*) realloc ( (my_mol1->chums), sizeof(int)*my_mol1->nchums );
my_mol1->chums[my_mol1->nchums-1]=chum2;

my_mol2->nchums++;
my_mol2->chums= (int*) realloc ( (my_mol2->chums), sizeof(int)*my_mol2->nchums );
my_mol2->chums[my_mol2->nchums-1]=chum1;

return FALSE;
}



/**Calculate the difference between two molecules
@param one molecule1
@param two molecule2
@return the difference
*/
double mol_diff(struct mol *one,struct mol *two)
{
int i;
double total=0.0;
struct vec delta;


	for (i=0;i<one->atoms;i++)
	{
		if (i>=two->atoms) return -1.0;
		cpy_vec(&delta,&(one->atom[i]));
		sub_vec(&delta,&(two->atom[i]));
		total+=mod_vec(&delta);
	}
return total;
}

/**Calculate the mass of a molecule
@param my_mol the molecule
*/
void mol_set_mass(struct mol *my_mol)
{
int i;
double total=0;
//		printf("%d |%s|\n",i,my_mol->element[0]);
//		getchar();

	for (i=0;i<my_mol->atoms;i++)
	{
		total+=ret_ram(my_mol->element[i]);
		if (i>63)
		{
			total+=2.0;
			if (i==(my_mol->atoms-1)) total+=1.0;
		}
	}
my_mol->mass=total;
}


///taken out of the main in need of a rewrite
void mol_get_ang(struct mol *my_mol,struct mol *my_mol2)
{
}

/**Calculate the radius of a molecule using its center of mass as the center of the circle.
*/
void mol_cal_r(struct mol *my_mol)
{
int i;
struct vec temp;
double max=0;
my_mol->r=0;
	for (i=0;i<(my_mol->atoms);i++)
	{
		cpy_vec(&(temp),&(my_mol->cog));
		sub_vec(&(temp),&(my_mol->atom[i]));
		//print_vec(&temp);
		max=mod_vec(&temp);
		if (max>my_mol->r) my_mol->r=max;
	}


}


///Rotate a molecule around an axis
void mol_rot(struct mol *mymol,struct vec *unit,struct vec *base, double ang)
{
	int ii;
	for (ii=0;ii<mymol->atoms;ii++)
	{
		rot_vec(&(mymol->atom[ii]),unit,base, ang);
	}
}

/**Add a number of atoms to a molecule.
@param my_mol the molecule on which to operate on 
@param atoms the number of atoms to add
*/
void mol_resize(struct mol *my_mol,int atoms)
{

	int ii;
	int oldlen=my_mol->atoms;
//printf("a\n");
//printf("b\n");

	my_mol->atom = (struct vec*) realloc(my_mol->atom,my_mol->atoms*sizeof(struct vec));
	my_mol->vel = (struct vec*) realloc(my_mol->vel,my_mol->atoms*sizeof(struct vec));
	my_mol->element = (char**) realloc(my_mol->element,my_mol->atoms*sizeof(char *));
//printf("c\n");

	for (ii=oldlen;ii<my_mol->atoms;ii++)
	{
		set_vec(&(my_mol->atom[ii]),-1,-1,-1);
		set_vec(&(my_mol->vel[ii]),-100,-100,-100);
		my_mol->element[ii] = (char*) malloc(4*sizeof(char ));

	}
//printf("d\n");


	my_mol->chainid = (char*) realloc(my_mol->chainid,sizeof(char)*my_mol->atoms);
	my_mol->dist = (char*) realloc(my_mol->dist,sizeof(char)*my_mol->atoms);
	my_mol->branch = (char*) realloc(my_mol->branch,sizeof(char)*my_mol->atoms);
//printf("e\n");
//my_mol->atoms--;
}

void mol_add_atom_to(struct mol *my_mol,char * name,struct vec *xyz)
{
int pos=my_mol->atoms;
mol_resize(my_mol,1);

strcpy(my_mol->element[pos],name);
my_mol->dist[pos]=' ';
my_mol->branch[pos]=' ';
my_mol->chainid[pos]=' ';

cpy_vec(&(my_mol->atom[pos]),xyz);
}

/**Initialize a molecule.
@param my_mol the molecule on which to operate on 
@param atoms the number of atoms that the molecule should have
*/
void mol_init(struct mol *my_mol,int atoms)
{
	my_mol->atoms=atoms;

	set_vec(&(my_mol->cog),-1,-1,-1);
	set_vec(&(my_mol->cog_n),-1,-1,-1);

	my_mol->atom = (struct vec*) malloc(my_mol->atoms*sizeof(struct vec));
	my_mol->vel = (struct vec*) malloc(my_mol->atoms*sizeof(struct vec));
	my_mol->element = (char**) malloc(my_mol->atoms*sizeof(char *));
	int ii;
	for (ii=0;ii<my_mol->atoms;ii++)
	{
		set_vec(&(my_mol->atom[ii]),-1,-1,-1);
		set_vec(&(my_mol->vel[ii]),-100,-100,-100);
		my_mol->element[ii] = (char*) malloc(4*sizeof(char ));

	}
	my_mol->resid = (char*) malloc(4*sizeof(char ));
	my_mol->chainid = (char*) malloc(sizeof(char)*my_mol->atoms);
	my_mol->dist = (char*) malloc(sizeof(char)*my_mol->atoms);
	my_mol->branch = (char*) malloc(sizeof(char)*my_mol->atoms);

	my_mol->chums=NULL;
	my_mol->nchums=0;
	my_mol->r=0;
	my_mol->del=FALSE;
//	set_vec(&(my_mol->step),0.0,0.0,0.0);

}

/**Free all memory associated with a molecule.
@param my_mol the molecule on which to operate on 
*/
void mol_free(struct mol *my_mol)
{

int i;
	for (i=0;i<my_mol->atoms;i++)
	{
		free(my_mol->element[i]);

	}

	free(my_mol->element);
	free(my_mol->atom);
	free(my_mol->vel);
	free(my_mol->chums);
	free(my_mol->resid);
	free(my_mol->chainid);
	free(my_mol->dist);
	free(my_mol->branch);
}

/**Rotate the molecule around the x-axis.
@param in the molecule on which to operate on 
@param ang the angle to rotate by 
*/
void mol_rotx(struct mol *in,double ang)
{
double a;
a=((2.0*pi)/360.0)*ang;
struct vec temp;
int i;

for (i=0;i<in->atoms;i++)
{
	temp.x=in->atom[i].x*1.0 + in->atom[i].y*0.0 +in->atom[i].z*0.0;

	temp.y=in->atom[i].x*0.0 + in->atom[i].y*cos(a) +in->atom[i].z*(-sin(a));
	temp.z=in->atom[i].x*0.0 + in->atom[i].y*sin(a) +in->atom[i].z*cos(a);
//printf("%lf %lf %lf %lf\n",in->atom[i].x*0.0 + in->atom[i].y*cos(a) +in->atom[i].z*(-sin(a)),in->atom[i].x , in->atom[i].y ,in->atom[i].z);
//getchar();
	in->atom[i].x=temp.x;
	in->atom[i].y=temp.y;
	in->atom[i].z=temp.z;

}

}

/**Rotate the molecule around the y-axis.
@param in the molecule on which to operate on 
@param ang the angle to rotate by 
*/
void mol_roty(struct mol *in,double ang)
{
double a;
a=((2.0*pi)/360.0)*ang;
struct vec temp;
int i;
for (i=0;i<in->atoms;i++)
{
	temp.x=in->atom[i].x*cos(a)    + in->atom[i].y*0.0 + in->atom[i].z*sin(a);
	temp.y=in->atom[i].x*0.0       + in->atom[i].y*1.0 + in->atom[i].z*0.0;
	temp.z=in->atom[i].x*(-sin(a)) + in->atom[i].y*0.0 + in->atom[i].z*cos(a);
	in->atom[i].x=temp.x;
	in->atom[i].y=temp.y;
	in->atom[i].z=temp.z;
}
}

/**Rotate the molecule around the z-axis.
@param in the molecule on which to operate on 
@param ang the angle to rotate by 
*/
void mol_rotz(struct mol *in,double ang)
{
double a;
a=((2.0*pi)/360.0)*ang;
struct vec temp;
int i;

for (i=0;i<in->atoms;i++)
{
	temp.x=in->atom[i].x*cos(a) + in->atom[i].y*(-sin(a)) +in->atom[i].z*0.0;
	temp.y=in->atom[i].x*sin(a) + in->atom[i].y*cos(a) +in->atom[i].z*0.0;
	temp.z=in->atom[i].x*0.0 + in->atom[i].y*0.0 +in->atom[i].z*1.0;
	in->atom[i].x=temp.x;
	in->atom[i].y=temp.y;
	in->atom[i].z=temp.z;
}

}


/**Calculate the cog of a molecule using atoms 0..n
@param my_mol the molecule on which to operate on 
@param n the number of atom sum over
*/
void mol_cal_cog_n(struct mol *my_mol,int n)
{
int i;
int max=n;
set_vec(&(my_mol->cog_n),0,0,0);
if (max>my_mol->atoms) max=my_mol->atoms;
	for (i=0;i<max;i++)
	{
		add_vec(&(my_mol->cog_n),&(my_mol->atom[i]));
	}

div_vec(&(my_mol->cog_n),(double)(max));
}


/**Calculate the cog of a molecule using all atoms
@param my_mol the molecule on which to operate on 
*/
void mol_cal_cog(struct mol *my_mol)
{
int i;
set_vec(&(my_mol->cog),0,0,0);
set_vec(&(my_mol->cog_n),0,0,0);
int count=0;
	for (i=0;i<(my_mol->atoms);i++)
	{
		add_vec(&(my_mol->cog),&(my_mol->atom[i]));
		if (count<60)
		{
			add_vec(&(my_mol->cog_n),&(my_mol->atom[i]));
			count++;
		}
	}

div_vec(&(my_mol->cog),(my_mol->atoms));
div_vec(&(my_mol->cog_n),(double)(count));
}


void mol_zero(struct mol *my_mol,int n)
{
struct vec delta;
mol_cal_cog(my_mol);
mol_cal_cog_n(my_mol,n);

cpy_vec(&delta,&(my_mol->cog_n));
mul_vec(&delta,-1);
mol_shift(my_mol,&delta);
mol_cal_cog(my_mol);
}

/**Shift the molecule by a vector
@param in the molecule on which to operate on 
@param delta the vector by which to shift the molecule
*/
void mol_shift(struct mol *in,struct vec *delta)
{
int i;
for (i=0;i<in->atoms;i++)
{
	in->atom[i].x+=delta->x;
	in->atom[i].y+=delta->y;
	in->atom[i].z+=delta->z;
}
	add_vec(&(in->cog),delta);
	add_vec(&(in->cog_n),delta);

}


/**Copy a molecule
@param my_mol1 destination molecule 
@param my_mol2 source molecule
*/
void mol_cpy(struct mol *my_mol1,struct mol *my_mol2)
{
	cpy_vec(&(my_mol1->cog),&(my_mol2->cog));
	cpy_vec(&(my_mol1->cog_n),&(my_mol2->cog_n));
	int i;
	
	for (i=0;i<my_mol1->atoms;i++)
	{
		//printf("%d %d\n",i,my_mol1->atoms);
		cpy_vec(&(my_mol1->atom[i]),&(my_mol2->atom[i]));
		strcpy(my_mol1->element[i],my_mol2->element[i]);
		my_mol1->chainid[i]=my_mol2->chainid[i];
		my_mol1->dist[i]=my_mol2->dist[i];
		my_mol1->branch[i]=my_mol2->branch[i];
	}
	//printf("end %d %d\n",my_mol1->nchums,my_mol2->nchums);

	my_mol1->nchums=my_mol2->nchums;
	my_mol1->chums= (int*) realloc ( (my_mol1->chums), sizeof(int)*my_mol1->nchums );

	for (i=0;i<my_mol2->nchums;i++)
	{
		my_mol1->chums[i]=my_mol2->chums[i];
	}
	//printf("here1\n");

	my_mol1->nchums=my_mol2->nchums;
	my_mol1->atoms=my_mol2->atoms;
	my_mol1->r=my_mol2->r;
//	cpy_vec(&(my_mol1->step),&(my_mol2->step));
	my_mol1->del=my_mol2->del;
	strcpy(my_mol1->resid,my_mol2->resid);

}
