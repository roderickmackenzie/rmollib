/** @file pdbgenbox.c
	@brief This program takes one molecule and places it in a bigger box as may times as you want/will fit.
	@details Command line arguments:\n
	@arg	 --out		Output file\n
	@arg	 --norot	no rotation\n
	@arg	 --bcb		on a bcb lattice\n
	@arg	 --orig		origonal input file\n
	@arg	 --dx		dx shift\n
	@arg	 --dy		dy shift\n
	@arg	 --shif		shift the whole box\n
	@arg	 --random	enable random placement\n
	@arg	 --mindist	minimum separation distance of two atoms\n
	@arg	 --needed	number of molecules needed\n
	@arg	 --load		file to load starting point from\n
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
//#include <mcheck.h> 
#include "lib/rmolib.h"

void help()
{
printf("--out\t Output file\n");
printf("\t--norot\t\t no rotation\n");
printf("\t--bcb\t\t on a bcb lattice\n");
printf("\t--orig\t\t origonal input file\n");
printf("\t--dx\t\t dx shift\n");
printf("\t--dy\t\t dy shift\n");
printf("\t--dz\t\t dz shift\n");
printf("\t--shift\t\t shift\n");
printf("\t--random\t random\n");
printf("\t--mindist\t mindist\n");
printf("\t--needed\t needed\n");
printf("\t--maxmass\t maxmass\n");
printf("\t--load\t\t load\n");

}


/**
Random box generation
*/
void gen_random_box(struct system *origsystem,struct system *in,struct box *mybox,char *name,int mindist,int needed,double shift)
{
progress_clear();
set_progress_multi_line();
set_porgress_max(50);
set_progress_colored();



int i=0;
double dist;
struct system temp;
//printf("\n\nattoms !!!! %d\n\n",origsystem.mols[0].atoms);

system_cal_cog(origsystem);

system_rescale(origsystem);
system_cal_r(origsystem);



system_init(&temp);

struct vec temp2;

int added=0;
int tried=0;

int ii;
int reject=FALSE;
int last=0;
int a;
in->size.x=mybox->xmax-mybox->xmin;
in->size.y=mybox->xmax-mybox->ymin;
in->size.z=mybox->xmax-mybox->zmin;

do
{
	struct vec delta;
	//cpy_mol(&temp,&origsystem.mols[0]);
//save_file("./mol.pdb",shift,&origsystem);
//printf("1\n");
//getchar();
	system_cpy(&temp,origsystem);

//save_file("./mol.pdb",shift,&origsystem);
//printf("2\n");
//getchar();

	delta.x=get_random_number(mybox->xmax);
	delta.y=get_random_number(mybox->ymax);
	delta.z=get_random_number(mybox->zmax);

	system_rotx(&temp,get_random_number(360));
	system_roty(&temp,get_random_number(360));
	system_rotz(&temp,get_random_number(360));
	system_shift(&temp,&delta);

	system_cal_cog(&temp);
	system_cal_r(&temp);


	struct vec test;

	for (i=0;i<temp.length;i++)
	{
		for (a=0;a<temp.mols[i].atoms;a++)	//loop over all the attoms in the mols
		{
			if (temp.mols[i].atom[a].x>mybox->xmax) reject=TRUE;
			if (temp.mols[i].atom[a].y>mybox->ymax) reject=TRUE;
			if (temp.mols[i].atom[a].z>mybox->zmax) reject=TRUE;

			if (temp.mols[i].atom[a].x<mybox->xmin) reject=TRUE;
			if (temp.mols[i].atom[a].y<mybox->ymin)	reject=TRUE;
			if (temp.mols[i].atom[a].z<mybox->zmin)	reject=TRUE;
		}

		if (reject==TRUE) break;

	}

	if (reject==FALSE)
	{
		//Now I have to decide if I want to reject the mol or not.
		for (i=0;i<temp.length;i++)
		{

			for (ii=0;ii<in->length;ii++)		//loop over mols
			{
				cpy_vec(&test,&(in->mols[ii].cog));
				sub_vec(&test,&(temp.mols[i].cog));
				dist=mod_vec(&test);
				if (dist<(temp.mols[i].r*2+mindist))
				{

					for (a=0;a<temp.mols[i].atoms;a++)	//loop over all the attoms in the mols
					{
						//struct vec umklapp;
						//set_vec(&umklapp,0.0,0.0,0.0);


						int b;
						for (b=0;b<in->mols[ii].atoms;b++)	//loop over all the attoms in the mols
						{

							cpy_vec(&temp2,&(in->mols[ii].atom[b]));
							//add_vec(&temp2,&umklapp);
							sub_vec(&temp2,&(temp.mols[i].atom[a]));
							dist=mod_vec(&temp2);

							if (dist<mindist)
							{
								//printf("Warning Reject\n");
								//getchar();
								reject=TRUE;
								break;
							}

						}

						if (reject==TRUE) break;

					}
				}
				if (reject==TRUE) break;

			}

			if (reject==TRUE) break;

		}
	}

	if (reject==FALSE)
	{
		struct vec null;
		set_vec(&null,0.0,0.0,0.0);
		add_to_system(in,&temp,&null);

		char text[100];
		sprintf(text,"added=%d needed=%d checked=%d effort=%d",added,needed,tried,tried-last);
		last=tried;

		set_progress_multi_line_text(text);
		text_progress((((double)added)/(double)needed));
		if ((added%100)==0) save_file(name,shift,in);
		added++;
	}
	tried++;
	reject=FALSE;


}while(added<needed);
text_progress_finish();




save_file(name,shift,in);

system_free(in);
printf("I have generated a total of %d molecules\n", added);
//printf("The total mass is %lf\n",origsystem.mols[0].mass*added);
}

/**
Standard cubic type lattice molecular generator.
*/
void gen_box(struct system *origsystem,struct system *in,struct box *mybox,struct vec *step, char *name, int norot,int bcb)
{
double new_dx_step=step->x;
if (bcb==TRUE) new_dx_step/=2.0;
double x;
double y;
double z;
int i=0;
int ytog=FALSE;
double ypush=0.0;
struct system temp;
system_init(&temp);
for (z=mybox->zmin+step->z/2.0;z<mybox->zmax-step->z/2.0;z+=step->z)
{

	for (x=mybox->xmin+step->x/2.0;x<mybox->xmax-step->x/2.0;x+=new_dx_step)
	{

		if (bcb==TRUE)
		{
			if (ytog==FALSE)
			{
				ypush=step->y/2.0;
				ytog=TRUE;
			}else
			{
				ypush=0.0;
				ytog=FALSE;
			}
		}

		for (y=mybox->ymin+step->y/2.0;y<mybox->ymax-step->y/2.0;y+=step->y)
		{	

		//printf("%lf %lf %lf\n",x,y,z);
			struct vec delta;
			set_vec(&delta,x,y+ypush,z);
			system_cpy(&temp,origsystem);

			system_shift(&temp,&delta);

			struct vec null;
			set_vec(&null,0.0,0.0,0.0);
	
			add_to_system(in,&temp,&null);

			i++;

			if (i>=in->length)
			{
				printf("I have run out of memory!!!!\n");
				exit(0);
			}
		
		}
	}
	if (bcb==TRUE)
	{
		if (ytog==FALSE)
		{
			ypush=step->y/2.0;
			ytog=TRUE;
		}else
		{
			ypush=0.0;
			ytog=FALSE;
		}
	}
}

struct vec min;
struct vec one;
set_vec(&one,-5.0,-5.0,-5.0);
system_get_min(in,&min);
add_vec(&min,&one);

if (min.x>0.0) min.x=0.0;
if (min.y>0.0) min.y=0.0;
if (min.z>0.0) min.z=0.0;

mul_vec(&min,-1.0);
system_shift(in,&min);




system_set_box(in);

in->size.x+=5.0;
in->size.y+=5.0;
in->size.z+=5.0;

save_file(name,0.0,in);
printf("%d",in->length);
system_free(in);
printf("I have generated a total of %d molecules\n", i);
}



int main(int argc, char* argv[])
{
rnd_init ();
int bcb=FALSE;
double shift=0.0;
struct system mainsystem;
struct system origsystem;
system_init(&mainsystem);
system_init(&origsystem);

	int norot=FALSE;

	if (scanarg( argv,argc,"--norot")) norot=TRUE;

	if (scanarg( argv,argc,"--bcb")) bcb=TRUE;

	if (get_arg_plusone_pos( argv,argc,"--orig")!=-1)
	{
		//double box_size=0.0;
		load_file(get_arg_plusone( argv,argc,"--orig"),&origsystem);
		system_post_load_fix_up(&origsystem);
		//printf("\n%d\n",orig[0].atoms);
		//getchar();
		//correct cog
		system_cal_cog(&origsystem);		
		//cal_cog(&origsystem.mols[0]);
		system_rescale(&origsystem);
		//shift_mol(&origsystem.mols[0],&mytemp);

		system_cal_cog(&origsystem);
		//cal_cog(&origsystem.mols[0]);
		//cal_r(&origsystem.mols[0]);
		struct vec null;
		set_vec(&null,0.0,0.0,0.0); 
		system_cal_r(&origsystem);
		
		//save_file("./fred.pdb",0,&origsystem);
		//cal_step(&origsystem.mols[0]);
		//get_ang(&origsystem.mols[0],0);
		//set_mass(&origsystem.mols[0]);
	}else
	{
		printf("No orig file.!\n");
		exit(0);
	}	

	struct box mybox;
	struct vec step;
	system_cal_step(&step,&origsystem);
	step.x+=2;
	step.y+=2;
	step.z+=2;

	printf("Box stepping at: ");
	print_vec(&step);

	if (scanarg( argv,argc,"--dx"))
	{
	step.x = atof(get_arg_plusone( argv,argc,"--dx"));
	}

	if (scanarg( argv,argc,"--dy"))
	{
	step.y = atof(get_arg_plusone( argv,argc,"--dy"));
	}

	if (scanarg( argv,argc,"--dz"))
	{
	step.z = atof(get_arg_plusone( argv,argc,"--dz"));
	}


	if (scanarg( argv,argc,"--shift"))
	{
		shift = atof(get_arg_plusone( argv,argc,"--shift"));
	}


	//For straight line molecules
	//mybox.xmin=0;
	//mybox.xmax=60;
	//mybox.ymin=0;
	//mybox.ymax=60;
	//mybox.zmin=0;
	//mybox.zmax=60;

	//For curlie molecules
	mybox.xmin=0;
	mybox.xmax=70;
	mybox.ymin=0;
	mybox.ymax=70;
	mybox.zmin=0;
	mybox.zmax=60;

	if (scanarg( argv,argc,"--xsize"))
	{

		mybox.xmax=atof(get_arg_plusone( argv,argc,"--xsize"));

	}else
	{
		printf("No --xsize");
		exit(0);
	}


	if (scanarg( argv,argc,"--ysize"))
	{


		mybox.ymax=atof(get_arg_plusone( argv,argc,"--ysize"));

	}else
	{
		printf("No --ysize");
		exit(0);
	}


	if (scanarg( argv,argc,"--zsize"))
	{

		mybox.zmax=atof(get_arg_plusone( argv,argc,"--zsize"));

	}else
	{
		printf("No --zsize");
		exit(0);
	}

	if (scanarg( argv,argc,"--random"))
	{
		int needed=100;
		int mindist=5;
		if (scanarg( argv,argc,"--mindist"))
		{
			mindist = atof(get_arg_plusone( argv,argc,"--mindist"));
		}else
		{
			printf("--mindist\n");
			exit(0);
		}

		if (scanarg( argv,argc,"--needed"))
		{
			needed = atof(get_arg_plusone( argv,argc,"--needed"));
		}else
		if (scanarg( argv,argc,"--maxmass"))
		{
			needed = atof(get_arg_plusone( argv,argc,"--maxmass"))/origsystem.mols[0].mass;
		}
		else
		{
			printf("no --needed or --maxmass\n");
			exit(0);
		}

		if (scanarg( argv,argc,"--load"))
		{
			load_file(get_arg_plusone( argv,argc,"--load"),&mainsystem);
			system_post_load_fix_up(&mainsystem);
		}
		else
		{		
			mainsystem.length=0;
			mainsystem.mols=NULL;
			//set_memory_mol(&mols,max_m,origsystem.mols[0].atoms);
		}

		gen_random_box(&origsystem,&mainsystem,&mybox, get_arg_plusone( argv,argc,"--out"),mindist,needed,shift);
	}else
	{	
		//system_free(&mainsystem);
		gen_box(&origsystem,&mainsystem,&mybox,&step, get_arg_plusone( argv,argc,"--out"),norot,bcb);
	}
rnd_free ();
return 0;

}
