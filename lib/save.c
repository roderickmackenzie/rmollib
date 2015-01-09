/** @file save.c
	@brief Routines to save pdb, pdbz, gro and com files
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include "minilzo.h"
#include "mol.h"

///Default call to save a system
void save_file(char * name,double shift,struct system *in)
{

char *type;
for (type=name+strlen(name);type>=name;type--)
{
	if (*type=='.')
	{
		type++;
		break;
	}
}


if (strcmp(type,"pdb")==0)
{
	dump_box_as_pdb(name,shift,in);
}else
if (strcmp(type,"gro")==0)
{
	dump_box_as_gro(in->mols,in->length, name ,in);
}else
if (strcmp(type,"pdbz")==0)
{
	save_pdb_comp(name,shift,in);
}else
{
	printf("I do not know what to do with this file\n");
	exit(0);
}

}

void left_align(char *out, int space,char *in)
{
	int len=strlen(in);
	int white=space-len;
	int i=0;
	strcat(out,in);
	for (i=0;i<white;i++)
	{
		strcat(out, " ");
	}

}


void save_pdb_frame(FILE *out,double shift,struct system *in)
{
struct expbuf buf;
expbuf_init(&buf);

build_pdb_frame(&buf,shift,in);
fwrite(buf.buf, buf.len, 1, out);
expbuf_free(&buf);
}



void save_pdb_frame_comp(FILE *out,double shift,struct system *in)
{
int i;
int r;
unsigned char* outbuf;
unsigned char* wrkmem;
lzo_uint out_len;

    if (lzo_init() != LZO_E_OK)
    {
        printf("internal error - lzo_init() failed !!!\n");
        printf("(this usually indicates a compiler bug - try recompiling\nwithout optimizations, and enable `-DLZO_DEBUG' for diagnostics)\n");
    }

unsigned char* buf=NULL;


int len=0;


for (i=0;i<in->length;i++)
{
	if (in->mols[i].del==FALSE)
	{
		len+=sizeof(float)*in->mols[i].atoms;	//x
		len+=sizeof(float)*in->mols[i].atoms;	//y
		len+=sizeof(float)*in->mols[i].atoms;	//z
		len+=sizeof(char)*in->mols[i].atoms;	//dist
		len+=sizeof(char)*in->mols[i].atoms;	//branch
		len+=sizeof(char)*in->mols[i].atoms;	//chain id
		len+=sizeof(char)*4*in->mols[i].atoms;	//element
		len+=sizeof(int);	//number of attoms
		len+=sizeof(char)*4;	//resid
		len+=sizeof(int);	//model
	}
}

len+=sizeof(int);		//how many mols are in the segment

buf = (unsigned char *) malloc(len);
outbuf = (unsigned char *) malloc(len);
build_z_frame(buf,shift,in);
//printf("comp\n");

wrkmem = (unsigned char *) malloc(LZO1X_1_MEM_COMPRESS);


r = lzo1x_1_compress((unsigned char* )buf,(lzo_uint)len,outbuf,&out_len,wrkmem);

if (r != LZO_E_OK)
	//printf("compressed %lu bytes into %lu bytes deflated by %lf\n",
   // (unsigned long) buf.len, (unsigned long) out_len, 100.0-((double)out_len/(double)buf.len)*100.0 );
{
/* this should NEVER happen */
printf("internal error - compression failed: %d\n", r);
}
//printf("write\n");
struct comp_header head;
head.len=(unsigned long) out_len;
head.origlen=(unsigned long) len;
head.xsize=in->size.x;
head.ysize=in->size.y;
head.zsize=in->size.z;
head.time=in->time;

fwrite(&head, sizeof(struct comp_header), 1, out);
//fwrite(buf, (unsigned long) len, 1, out);
fwrite(outbuf, (unsigned long) out_len, 1, out);

//printf("thought I wrote%ld\n",len);
free(buf);
free(wrkmem);
free(outbuf);
}


///Build compressed pdb frame in binary format
void build_z_frame(unsigned char *out,double shift,struct system *in)
{

int i;
int ii;
unsigned char* pos=out;
int mols=0;

for (i=0;i<in->length;i++)
{

	if (in->mols[i].del==FALSE)
	{
		mols+=1;
	}
}

memcpy(pos, &(mols), sizeof(int));
pos+=sizeof(int);

for (i=0;i<in->length;i++)
{

	if (in->mols[i].del==FALSE)
	{

		memcpy(pos, &(in->mols[i].atoms), sizeof(int));
		pos+=sizeof(int);

		memcpy(pos, &(in->model), sizeof(int));
		pos+=sizeof(int);

		memcpy(pos, in->mols[i].resid, sizeof(char)*4);
		pos+=sizeof(char)*4;


		float f;
		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			f=in->mols[i].atom[ii].x;
			memcpy(pos,&f,sizeof(float));
			pos+=sizeof(float);
		}

		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			f=in->mols[i].atom[ii].y;
			memcpy(pos,&f,sizeof(float));
			pos+=sizeof(float);
		}


		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			f=in->mols[i].atom[ii].z;
			memcpy(pos,&f,sizeof(float));
			pos+=sizeof(float);
		}


		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(pos,&(in->mols[i].dist[ii]),sizeof(char));
			pos+=sizeof(char);
		}


		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(pos,&(in->mols[i].branch[ii]),sizeof(char));
			pos+=sizeof(char);
		}


		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(pos,&(in->mols[i].chainid[ii]),sizeof(char));
			pos+=sizeof(char);
		}

		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(pos,in->mols[i].element[ii],4);
			pos+=sizeof(char)*4;
		}



	}
}
//			pos-=sizeof(char)*3;
//			*(pos++)='r';
//			*(pos++)='o';
//			*(pos++)='d';
//printf("written %d\n",pos-out);
}


///build a normal pdb file.
void build_pdb_frame(struct expbuf *out,double shift,struct system *in)
{
char line[1000];
int num;
num=1;
int mol;
mol=1;

if (in->time!=-1.0)
{
	sprintf(line,"TITLE     Protein t=   %lf\n",in->time);
	expbuf_add(out,line);
}

if (in->size.x!=-1) sprintf(line,"CRYST1  %7.3f  %7.3f  %7.3f  90.00  90.00  90.00 P 1           1\n",in->size.x,in->size.y,in->size.z);
expbuf_add(out,line);

sprintf(line,"MODEL        %d\n",in->model);
expbuf_add(out,line);

int i;
int ii;
char build[500];
char pretected_resid[100];
//printf("the length is now %d\n",in->length);
for (i=0;i<in->length;i++)
{
				//printf("here %d %d %d\n",mol,in->length,in->mols[i].atoms);
				//getchar();
	//printf("%d\n",i);

	if (in->mols[i].del==FALSE)
	{

		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			strcpy(pretected_resid,in->mols[i].resid);
			pretected_resid[3]=0;		//terminate the string at 3 whatever
			//if (in->mols[i].branch[ii]==0) in->mols[i].branch[ii]=' ';

			sprintf(build,"%-6s%5d%c%2s%c%c%c%s%c%c%4d%c   %8.3lf%8.3lf%8.3lf  1.00  0.00\n","ATOM",num,' ',in->mols[i].element[ii],in->mols[i].dist[ii],in->mols[i].branch[ii],' ',pretected_resid,' ',in->mols[i].chainid[ii],mol,' ',in->mols[i].atom[ii].x+shift,in->mols[i].atom[ii].y+shift,in->mols[i].atom[ii].z+shift);
			//sprintf(line,"%s\n",build);
			expbuf_add(out,build);

			if (num>=99999)
			{
				//printf("Warning pdb file over 99999 entries!\n");
				num=0;
			}else
			{
				num++;
			}
		}

		if (mol>=9999)
		{
			//printf("Warning pdb file over 99999 entries!\n");
			mol=0;
		}else
		{
			mol++;
		}
	}
}


sprintf(line,"TER\nENDMDL\n");
expbuf_add(out,line);

}


///Dump a system to a file in a pdb format.
void dump_box_as_pdb(char *name,double shift,struct system *in)
{
FILE *out;
out=fopen(name,"wb");//fopen(name,"w");

save_pdb_frame(out,shift,in);

fclose(out);


}


void save_pdb_comp(char *name,double shift,struct system *in)
{
FILE *out;
out=fopen(name,"wb");//fopen(name,"w");

save_pdb_frame_comp(out,shift,in);

fclose(out);


}


void dump_box_as_gro(struct mol *myarray,int length, char *name,struct system *mysystem)
{
int i;
int ii;
int num;
num=1;
int mol;
mol=1;
FILE *out;
int tot=0;
for (i=0;i<length;i++)
{
	if (myarray[i].del==FALSE) tot+=myarray[i].atoms;
}

out=fopen(name,"w");
fprintf(out,"Generated by rods gro generator\n");
fprintf(out,"%5d\n",tot);



for (i=0;i<length;i++)
{
//printf("%d\n",i);
	if (myarray[i].del==FALSE)
	{

		for (ii=0;ii<myarray[i].atoms;ii++)
		{
			char asem[10];
			sprintf(asem,"%3s%c%c",myarray[i].element[ii],myarray[i].dist[ii],myarray[i].branch[ii]);
			remove_white_space(asem);
			if (myarray[i].vel[ii].z!=-100.0)
			{

				fprintf(out,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",mol,myarray[i].resid,asem,num,myarray[i].atom[ii].x/10.0,myarray[i].atom[ii].y/10.0,myarray[i].atom[ii].z/10.0,myarray[i].vel[ii].x,myarray[i].vel[ii].y,myarray[i].vel[ii].z);

			}else
			{
				fprintf(out,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",mol,myarray[i].resid,asem,num,myarray[i].atom[ii].x/10.0,myarray[i].atom[ii].y/10.0,myarray[i].atom[ii].z/10.0);
			}



			num++;
			if (num>99999) num=0;
		}

		mol++;
		if (mol>99999) mol=0;
	}

}

fprintf(out,"%10.5f%10.5f%10.5f\n",mysystem->size.x/10.0,mysystem->size.y/10.0,mysystem->size.z/10.0);
fclose(out);
}

///Dump a com file
void dump_com_file(struct mol *my_vec1,struct mol *my_vec2, char *name,struct vec *cog,char *gaus_base_file)
{
}


///Dump a pair of molecules as a normal pdb file
void dump_pair(struct mol *my_vec1,struct mol *my_vec2,char *name,struct vec *shift)
{

int num;
num=1;
int mol;
mol=1;
FILE *out;
out=fopen(name,"w");
fprintf(out,"MODEL        1\n");
int ii;
char build[1000];


if (my_vec1!=0)
{
//printf("one\n");
	for (ii=0;ii<my_vec1->atoms;ii++)
	{
sprintf(build,"%-6s%5d%c%2s%c%c%c%3s%c%c%4d%c   %8.3lf%8.3lf%8.3lf  1.00  0.00","ATOM",num,' ',my_vec1->element[ii],my_vec1->dist[ii],my_vec1->branch[ii],' ',my_vec1->resid,' ',my_vec1->chainid[ii],mol,' ',my_vec1->atom[ii].x+shift->x,my_vec1->atom[ii].y+shift->y,my_vec1->atom[ii].z+shift->z);
		fprintf(out,"%s\n",build);

		num++;
	}
	mol++;
}

if (my_vec2!=0)
{


	for (ii=0;ii<my_vec2->atoms;ii++)
	{
sprintf(build,"%-6s%5d%c%2s%c%c%c%3s%c%c%4d%c   %8.3lf%8.3lf%8.3lf  1.00  0.00","ATOM",num,' ',my_vec2->element[ii],my_vec2->dist[ii],my_vec2->branch[ii],' ',my_vec2->resid,' ',my_vec2->chainid[ii],mol,' ',my_vec2->atom[ii].x+shift->x,my_vec2->atom[ii].y+shift->y,my_vec2->atom[ii].z+shift->z);
		fprintf(out,"%s\n",build);

		num++;
	}
}

fprintf(out,"TER\nENDMDL\n");
fclose(out);
}

