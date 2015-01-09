/** @file system.c
	@brief Algorithms to manipulate entire molecular systems.
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
#include "molops.h"
#include "progress.h"



/**Count the number of a particular number of residuals
@param in the system to process
@param name name of the residue
@return the number of residuals
*/
int system_count_resid(struct system * in,char *name)
{
int i;
int count=0;
for (i=0;i<in->length;i++)
{
	if (strcmp(in->mols[i].resid,name)==0)  count++;
}
return count;
}

///Calculate the cog of a system using atoms 0..n 
void system_cal_cog_n(struct system *in,int n)
{
int i;
for (i=0;i<in->length;i++)
{
mol_cal_cog_n(&(in->mols[i]),n);
}

}
void system_find_resids(struct system * in)
{
int *list;
list = NULL;

int *list_count;
list_count = NULL;

int len=0;

int i;
int ii;
int found=FALSE;
for (i=0;i<in->length;i++)
{

	found=FALSE;

	for (ii=0;ii<len;ii++)
	{
		if (strcmp(in->mols[i].resid,in->mols[list[ii]].resid)==0)
		{
			found=TRUE;
			list_count[ii]++;
		}
	}

	if (found==FALSE)
	{
		list=(int *)realloc(list,sizeof(int)*(len+1));
		list_count=(int *)realloc(list_count,sizeof(int)*(len+1));
		list_count[len]=1;
		list[len]=i;
		len++;
	}	
}
for (i=0;i<len;i++)
{
	printf("%s %d\n",in->mols[list[i]].resid,list_count[i]);
}
free(list);
free(list_count);
}


//the point of this function is to calculate angles COGs etc just after the load.
void system_post_load_fix_up(struct system *in)
{
int i=0;

for (i=0;i<in->length;i++)
{
	mol_cal_cog(&(in->mols[i]));
	mol_cal_r(&(in->mols[i]));
	mol_get_ang(&(in->mols[i]),0);
	mol_set_mass(&(in->mols[i]));
}

}

/**Assigns a group number to residules, this is used primeraly for polymers of a set length.
@param mysystem the system to process
@param len length of groups to make
*/
void system_set_group(struct system *mysystem,int len)
{
int i;
int count=0;
int group=0;
for (i=0;i<mysystem->length;i++)
{
	//if (strcmp(mysystem->mols[i].resid,resid)==0)
	//{
		mysystem->mols[i].group=group;
		//printf("Just added to group %d\n",group);
		count++;
		//printf("Coiunt at end of add is %d\n",count);

		if (count>=len)
		{
			count=0;
			group++;
		}



		//getchar();
	//}
}

}

void rot_system(struct system *mysystem,struct vec *unit,struct vec *base, double ang)
{
	int ii;
	for (ii=0;ii<mysystem->length;ii++)
	{
		mol_rot(&(mysystem->mols[ii]),unit,base, ang);
	}
}


double system_cal_density(struct system *in)
{
int i;
int ii;

text_progress_start("Calculating density....");
set_porgress_max(40);
set_progress_colored();


double x_min=in->systembox.xmin+in->size.x*0.2;
double x_max=in->systembox.xmax-in->size.x*0.2;

double x_len=x_max-x_min;

double y_min=in->systembox.ymin+in->size.y*0.2;
double y_max=in->systembox.ymax-in->size.y*0.2;

double y_len=y_max-y_min;

double z_min=in->systembox.zmin+in->size.z*0.2;
double z_max=in->systembox.zmax-in->size.z*0.2;

double z_len=z_max-z_min;
double mr=1.660538782e-24;	//mass

double total=0;

for (i=0;i<in->length;i++)
{
if (i%100==0)text_progress((((double)i)/(double)in->length));



	for (ii=0;ii<in->mols[i].atoms;ii++)
	{
		if ((in->mols[i].atom[ii].x>=x_min)&&(in->mols[i].atom[ii].x<=x_max))
		{
			if ((in->mols[i].atom[ii].y>=y_min)&&(in->mols[i].atom[ii].y<=y_max))
			{
				if ((in->mols[i].atom[ii].z>=z_min)&&(in->mols[i].atom[ii].z<=z_max))
				{
					total+=ret_ram(in->mols[i].element[ii]);
				}
			}
		}
	
	}

}

text_progress_finish();
return total*mr/(y_len*x_len*z_len*1e-8*1e-8*1e-8);
}

///Shift the system by a set vector
void system_shift(struct system *in,struct vec *delta)
{
int i;
	for (i=0;i<in->length;i++)
	{
		mol_shift(&(in->mols[i]),delta);
	}
}

void system_rotx(struct system *in,double ang)
{
int i;
	for (i=0;i<in->length;i++)
	{
		mol_rotx(&(in->mols[i]),ang);
	}
}

///Rotate system on y axis
void system_roty(struct system *in,double ang)
{
int i;
	for (i=0;i<in->length;i++)
	{
		mol_roty(&(in->mols[i]),ang);
	}
}

///Rotate system on z axis
void system_rotz(struct system *in,double ang)
{
int i;
	for (i=0;i<in->length;i++)
	{
		mol_rotz(&(in->mols[i]),ang);
	}
}	

///Truncate the system velocities.  Supply a number you want it to round to i.e. 1.0 or 0.0001.
void system_truncate_v(struct system *in,double v)
{
int i;
int ii;

	for (i=0;i<in->length;i++)
	{
		for (ii=0;ii<(in->mols[i].atoms);ii++)
		{
			in->mols[i].vel[ii].x=(double)((int)(in->mols[i].vel[ii].x/v))*v;
			in->mols[i].vel[ii].y=(double)((int)(in->mols[i].vel[ii].y/v))*v;
			in->mols[i].vel[ii].z=(double)((int)(in->mols[i].vel[ii].z/v))*v;

		}
	}
}


void system_get_min(struct system *in,struct vec *ret)
{
int i;
int ii;
double x=in->mols[0].atom[0].x;
double y=in->mols[0].atom[0].y;
double z=in->mols[0].atom[0].z;

	for (i=0;i<in->length;i++)
	{
		for (ii=0;ii<(in->mols[i].atoms);ii++)
		{
			if (in->mols[i].atom[ii].x<x) x=in->mols[i].atom[ii].x;
			if (in->mols[i].atom[ii].y<y) y=in->mols[i].atom[ii].y;
			if (in->mols[i].atom[ii].z<z) z=in->mols[i].atom[ii].z;

		}
	}

set_vec(ret,x,y,z);
}

///Compares two systems and gives an integer describing the difference between the two systems.
///if resid is set to "" the whole system is compared.
///The function assumes that the molecules are placed in the same order in the file.
double system_diff(struct system *one,struct system *two,char *resid)
{
double tot=0.0;
double ret=0.0;
int i;
int compare;
int max=one->length;
if (two->length<max) max=two->length;

	for (i=0;i<max;i++)
	{
		compare=FALSE;

		if (strcmp(one->mols[i].resid,resid)==0)
		if (strcmp(two->mols[i].resid,resid)==0)
		{
			compare=TRUE;
		}

		if (strcmp(resid,"")==0) compare=TRUE;

		if (compare==TRUE)
		{
			ret=mol_diff(&(one->mols[i]),&(two->mols[i]));
			if (ret<0) return -1;
			tot+=ret;
		}
	}
return tot;
}

void system_rescale(struct system *in)
{
int i;
int ii;

	for (i=0;i<in->length;i++)
	{
		for (ii=0;ii<(in->mols[i].atoms);ii++)
		{
			sub_vec(&(in->mols[i].atom[ii]),&(in->cog));
		}
	}


}

void system_cal_cog(struct system *in)
{
int i;
int ii;
struct vec temp;

set_vec(&temp,0.0,0.0,0.0);

int total=0;
	for (i=0;i<in->length;i++)
	{
		for (ii=0;ii<(in->mols[i].atoms);ii++)
		{
			add_vec(&(temp),&(in->mols[i].atom[ii]));
		}
		mol_cal_cog(&(in->mols[i]));
		total+=in->mols[i].atoms;
	}

div_vec(&temp,total);

cpy_vec(&(in->cog),&temp);
	
};

void system_cal_r(struct system *in)
{
int i;
int ii;
struct vec temp;


double max=0;
in->r=0;

	for (i=0;i<in->length;i++)
	{
		for (ii=0;ii<(in->mols[i].atoms);ii++)
		{
			cpy_vec(&(temp),&(in->cog));
			sub_vec(&(temp),&(in->mols[i].atom[ii]));
			max=mod_vec(&temp);
			if (max>in->mols[i].r) in->mols[i].r=max;
		}

		mol_cal_r(&(in->mols[i]));

	}


};

void system_cal_step(struct vec *step,struct system *in)
{
double min_x=in->mols[0].atom[0].x;
double max_x=in->mols[0].atom[0].x;
double min_y=in->mols[0].atom[0].y;
double max_y=in->mols[0].atom[0].y;
double min_z=in->mols[0].atom[0].z;
double max_z=in->mols[0].atom[0].z;

int i;
int ii;
	for (i=0;i<in->length;i++)
	{
		for (ii=0;ii<(in->mols[i].atoms);ii++)
		{
			if (in->mols[i].atom[ii].x<min_x) min_x=in->mols[i].atom[ii].x;
			if (in->mols[i].atom[ii].x>max_x) max_x=in->mols[i].atom[ii].x;

			if (in->mols[i].atom[ii].y<min_y) min_y=in->mols[i].atom[ii].y;
			if (in->mols[i].atom[ii].y>max_y) max_y=in->mols[i].atom[ii].y;

			if (in->mols[i].atom[ii].z<min_z) min_z=in->mols[i].atom[ii].z;
			if (in->mols[i].atom[ii].z>max_z) max_z=in->mols[i].atom[ii].z;
		}

	}

step->x=max_x-min_x;
step->y=max_y-min_y;
step->z=max_z-min_z;


}

/// Initialize the system.
void system_init(struct system* in)
{
in->model=0;
in->length=0;
in->mols=NULL;
set_vec(&(in->size),-1,-1,-1);
}

///Free all the memory given to the system.
void system_free(struct system *in)
{
int i;
for (i=0;i<in->length;i++)
{
	mol_free(&(in->mols[i]));	
}

free(in->mols);
system_init(in);
}

///Make a copy of a system
void system_cpy(struct system *out,struct system *in)
{
struct vec shift;

set_vec(&shift,0.0,0.0,0.0);

system_free(out);

out->r=in->r;

cpy_vec(&(out->cog),&(in->cog));

add_to_system(out,in,&shift);

}

struct mol * mol_expand_memory(struct system * in,int atoms)
{
in->length++;
in->mols = (struct mol *) realloc(in->mols, in->length*sizeof(struct mol));
mol_init(&(in->mols[in->length-1]),atoms);
return &(in->mols[in->length-1]);
}




void system_insert_mol(struct system *out,struct mol *in,struct vec *shift)
{
	mol_expand_memory(out,in->atoms);
	
	mol_cpy(&out->mols[out->length-1],in);
	struct vec temp;
	cpy_vec(&temp,&out->mols[out->length-1].cog_n);
	mul_vec(&temp,-1.0);
	mol_shift(&out->mols[out->length-1],&temp);
	mol_shift(&out->mols[out->length-1],shift);
	out->mols[out->length-1].del=FALSE;
}

///Add one system to another system
void add_to_system(struct system *out,struct system *in,struct vec *shift)
{
int i;

for (i=0;i<in->length;i++)
{
	mol_expand_memory(out,in->mols[i].atoms);
	mol_cpy(&out->mols[out->length-1],&in->mols[i]);
	mol_shift(&out->mols[out->length-1],shift);

	//printf("a rod %d %d\n",out->length,out->mols[out->length-1].atoms);

	//getchar();
}

}

///Delete molecules from the system with the name equal to resid
void system_del_mols(struct system *in,char *resid)
{
int i;

for (i=0;i<in->length;i++)
{

	if (strcmp(in->mols[i].resid,resid)==0) in->mols[i].del=TRUE;

}


}

///Delete molecules from the system with the name not equal to resid
void system_inv_del_mols(struct system *in,char *resid)
{
int i;

for (i=0;i<in->length;i++)
{

	if (strcmp(in->mols[i].resid,resid)!=0) in->mols[i].del=TRUE;

}


}


///Updated function try to make other two work
void system_set_box_now(struct system *in)
{
int i;
int ii;

double x_max=in->mols[0].atom[0].x;
double y_max=in->mols[0].atom[0].y;
double z_max=in->mols[0].atom[0].z;

double x_min=in->mols[0].atom[0].x;
double y_min=in->mols[0].atom[0].y;
double z_min=in->mols[0].atom[0].z;

for (i=0;i<in->length;i++)
{
	for (ii=0;ii<in->mols[i].atoms;ii++)
	{
		if (in->mols[i].del==FALSE)
		{
			if (in->mols[i].atom[ii].x>x_max) x_max=in->mols[i].atom[ii].x;
			if (in->mols[i].atom[ii].y>y_max) y_max=in->mols[i].atom[ii].y;
			if (in->mols[i].atom[ii].z>z_max) z_max=in->mols[i].atom[ii].z;

			if (in->mols[i].atom[ii].x<x_min) x_min=in->mols[i].atom[ii].x;
			if (in->mols[i].atom[ii].y<y_min) y_min=in->mols[i].atom[ii].y;
			if (in->mols[i].atom[ii].z<z_min) z_min=in->mols[i].atom[ii].z;
		}
	}
	//printf("%lf %lf %lf %lf %lf %lf\n",x_min,y_min,z_min,x_max,y_max,z_max);
}

in->systembox.xmin=0.0;//x_min;
in->systembox.xmax=0.0;//x_max;
in->systembox.ymin=0.0;//y_min;
in->systembox.ymax=0.0;//y_max;
in->systembox.zmin=0.0;//z_min;
in->systembox.zmax=0.0;//z_max;


in->size.x=x_max-x_min;
in->size.y=y_max-y_min;
in->size.z=z_max-z_min;
struct vec min;
set_vec(&min,x_min,y_min,z_min);
mul_vec(&min,-1.0);
system_shift(in,&min);

}

void system_set_box(struct system *in)
{
int i;
int ii;

double x_max=in->mols[0].atom[0].x;
double y_max=in->mols[0].atom[0].y;
double z_max=in->mols[0].atom[0].z;

double x_min=in->mols[0].atom[0].x;
double y_min=in->mols[0].atom[0].y;
double z_min=in->mols[0].atom[0].z;

for (i=0;i<in->length;i++)
{
	for (ii=0;ii<in->mols[i].atoms;ii++)
	{
		if (in->mols[i].del==FALSE)
		{
			if (in->mols[i].atom[ii].x>x_max) x_max=in->mols[i].atom[ii].x;
			if (in->mols[i].atom[ii].y>y_max) y_max=in->mols[i].atom[ii].y;
			if (in->mols[i].atom[ii].z>z_max) z_max=in->mols[i].atom[ii].z;

			if (in->mols[i].atom[ii].x<x_min) x_min=in->mols[i].atom[ii].x;
			if (in->mols[i].atom[ii].y<y_min) y_min=in->mols[i].atom[ii].y;
			if (in->mols[i].atom[ii].z<z_min) z_min=in->mols[i].atom[ii].z;
		}
	}
	//printf("%lf %lf %lf %lf %lf %lf\n",x_min,y_min,z_min,x_max,y_max,z_max);
}

in->systembox.xmin=x_min;
in->systembox.xmax=x_max;
in->systembox.ymin=y_min;
in->systembox.ymax=y_max;
in->systembox.zmin=z_min;
in->systembox.zmax=z_max;


in->size.x=x_max;
in->size.y=y_max;
in->size.z=z_max;


}





///This will take one system with periodic boundary conditions and expand it to make a bigger block from it.
///If tidy is set it will delete atoms with a center out side of the box.
void system_expand(struct system *out,struct system *in,double size_x,double size_y,double size_z, int tidy)
{
//int x;
//int y;
//int z;
int i;
//int ii;

//progress_clear();
out->time=in->time;
//text_progress_start("Working");

set_porgress_max(40);
set_progress_colored();
double step_x=in->size.x;
double step_y=in->size.y;
double step_z=in->size.z;
step_x=round(step_x);
step_y=round(step_y);
step_z=round(step_z);

//int max_x=(int)(size_x/step_x)+1;
//int max_y=(int)(size_y/step_y)+1;
//int max_z=(int)(size_z/step_z)+1;

//printf("%d %d %d\n",max_x,max_y,max_z);

int count=0;
//int max=max_x*max_y*max_z;
double xpos=0;
double ypos=0;
double zpos=0;
int xloop=0;
int yloop=0;
int zloop=0;
do
{

	yloop=0;
	do
	{

		zloop=0;
		do
		{	

			struct vec delta;
			set_vec(&delta,xpos,ypos,zpos);
			//printf("a %.10lf %.10lf %.10lf\n",((double)x)*step_x,((double)y)*step_y,((double)z)*step_z);
			//printf("b %.10lf %.10lf %.10lf\n",((double)x)*in->size.x,((double)y)*in->size.y,((double)z)*in->size.z);
			add_to_system(out,in,&delta);
			//text_progress((((double)count)/(double)(size_z/xpos)*(size_z/xpos)*(size_z/xpos)));
			count++;

		zpos+=step_z;
		zloop++;
		}while(zpos<(size_z));
		zpos=0.0;

	ypos+=step_y;
	yloop++;
	}while(ypos<(size_y));
	ypos=0.0;

xpos+=step_x;
xloop++;
}while(xpos<(size_x));

//system_set_box_funct(out);
//system_set_box(out);
//printf("fin\n");

//double x_max=out->funcbox.xmin+size_x;
//double y_max=out->funcbox.ymin+size_y;
//double z_max=out->funcbox.zmin+size_z;
//in->size.x=size_x;
//in->size.y=size_y;
//in->size.z=size_z;
//printf("%lf %lf %lf\n",x_max,y_max,z_max);
//printf("%lf %lf %lf\n",size_x,size_y,size_z);
//printf("%lf %lf %lf\n",in->systembox.xmin,in->systembox.ymin,in->systembox.zmin);

//getchar();
if (tidy==TRUE)
{
system_cal_cog_n(out,19);
//printf("Tiday\n");
for (i=0;i<out->length;i++)
{
//printf("%lf\n",out->mols[i].cog_n.x);
//getchar();
		if ((out->mols[i].cog_n.x>size_x)||(out->mols[i].cog_n.y>size_y)||(out->mols[i].cog_n.z>size_z))
		{
			//printf("Deleted\n");
			out->mols[i].del=TRUE;
		}

		if ((out->mols[i].cog_n.x<0.0)||(out->mols[i].cog_n.y<0.0)||(out->mols[i].cog_n.z<0.0))
		{
			//printf("Deleted\n");
			out->mols[i].del=TRUE;
		}

}
}
out->size.x=((double)xloop)*step_x;
out->size.y=((double)yloop)*step_y;
out->size.z=((double)zloop)*step_z;
//printf("%ld %ld %ld\n",xloop,yloop,zloop);
//system_set_box_funct(out);
//system_set_box(out);


//text_progress_finish();
}
