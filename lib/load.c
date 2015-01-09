/** @file load.c
	@brief Routines to load pdb, pdbz, gro and com files
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
#include "load.h"

/**Generic load function it will load pdb file, gro files and pdbz files.
*/
void load_file(char * name,struct system *in)
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
	//printf("rod= %d\n",in->length);
	load_pdb(name,in);

}else
if (strcmp(type,"gro")==0)
{
	in->length=0;
	load_gro(name,in);
}else
if (strcmp(type,"pdbz")==0)
{
	in->length=0;
	load_pdb_comp(name,in);
}else
{
	printf("I do not know what to do with this file\n");
	exit(0);
}

}

///Wrapper for bzero.
void set_null(char *in,int len)
{
bzero(in, len);
}





///I have no idea what this does.
int util_find_end_of(char *in,char *check)
{
int i;
int ii;
int found;
found=FALSE;
int pos=-1;
	for (i=0;i<strlen(in);i++)
	{
		int max=i+strlen(check);
		if ((max)>i+strlen(in))
		{
		max=i+strlen(in);
		}
	
		found=TRUE;
		for (ii=0;ii<strlen(check);ii++)
		{
			//printf("%c %c %d %d\n",in[i+ii],check[ii],i,ii);
			if (in[i+ii]!=check[ii])
			{
				found=FALSE;
				break;
			}
		}
		
		if (found==TRUE)
		{
			pos=i;
			break;
		}
}

return pos;
}

///Loads a normal uncompressed pdb frame from a file.  Returns FALSE is the read has not been successful.
int load_pdb_frame(FILE *file,struct system *in,int justgetinfo)
{
char buffer[1000];
struct expbuf buf;
expbuf_init(&buf);
int quit=FALSE;
do
{
	fgets(buffer, 1000, file);
	expbuf_add(&buf,buffer);

	if (mycmp(buffer,"END"))
	{
		quit=TRUE;
	}

	if (quit==FALSE)
	{
		if (feof(file)!=0)
		{
			return FALSE;
		}
	}
}while (quit==FALSE);
decode_pdb_frame(&buf,in,justgetinfo);
expbuf_free(&buf);

return TRUE;
}


///Load a compressed pdb file frame in the pdbz format
void load_pdb_frame_comp(FILE *file,struct system *in,int justgetinfo)
{
int r;
unsigned char* inbuf;
unsigned char* wrkmem;
lzo_uint out_len;

struct comp_header head;
unsigned char* outbuf;

int myread=fread(&head, sizeof (struct comp_header ), 1, file);

if (myread<=0) return;



if (justgetinfo==TRUE)
{
in->size.x=head.xsize;
in->size.y=head.ysize;
in->size.z=head.zsize;
in->time=head.time;
return;
}


wrkmem = (unsigned char *) malloc(LZO1X_1_MEM_COMPRESS);
outbuf = (unsigned char *) malloc((head.origlen+100)*sizeof(char));


inbuf = (unsigned char *) malloc(head.len*sizeof(unsigned char));
    if (lzo_init() != LZO_E_OK)
    {
        printf("internal error - lzo_init() failed !!!\n");
        printf("(this usually indicates a compiler bug - try recompiling\nwithout optimizations, and enable `-DLZO_DEBUG' for diagnostics)\n");
    }

fread(inbuf, sizeof(unsigned char)*head.len, 1, file);
    r = lzo1x_decompress(inbuf,head.len,(unsigned char *)outbuf,&out_len,NULL);
    if (r == LZO_E_OK && head.origlen == out_len)
    {
       // printf("inflateing %lf\n",
        //    100.0*((double) out_len/(double) head.len));
    }else
    {
        /* this should NEVER happen */
        printf("internal error - decompression failed: %d\n", r);
	exit(0);
    }


decode_z_frame(outbuf,in,justgetinfo);

in->size.x=head.xsize;
in->size.y=head.ysize;
in->size.z=head.zsize;
in->time=head.time;
//printf("%lf",in->size.x);
//getchar();

free(outbuf);
free(wrkmem);
free(inbuf);
}


///Decode a compressed pdbz frame.
void decode_z_frame(unsigned char *buf,struct system *in,double time)
{
int i;
int ii;
unsigned char* pos=buf;
int len;
system_free(in);
memcpy(&(len),pos, sizeof(int));
pos+=sizeof(int);

for (i=0;i<len;i++)
{
		int atoms=0;
		memcpy(&(atoms),pos, sizeof(int));
		pos+=sizeof(int);
		mol_expand_memory(in,atoms);

		memcpy(&(in->model),pos, sizeof(int));
		pos+=sizeof(int);

		memcpy(in->mols[i].resid, pos, sizeof(char)*4);
		pos+=sizeof(char)*4;


		float f;


		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(&f,pos,sizeof(float));
			in->mols[i].atom[ii].x=(double)f;
			pos+=sizeof(float);
		}

		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(&f,pos,sizeof(float));
			in->mols[i].atom[ii].y=(double)f;
			pos+=sizeof(float);
		}



		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(&f,pos,sizeof(float));
			in->mols[i].atom[ii].z=(double)f;
			pos+=sizeof(float);
		}




		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(&(in->mols[i].dist[ii]),pos,sizeof(char));
			pos+=sizeof(char);
		}


		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(&(in->mols[i].branch[ii]),pos,sizeof(char));
			pos+=sizeof(char);
		}


		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(&(in->mols[i].chainid[ii]),pos,sizeof(char));
			pos+=sizeof(char);
		}

		for (ii=0;ii<in->mols[i].atoms;ii++)
		{
			memcpy(in->mols[i].element[ii],pos,4);
			pos+=sizeof(char)*4;
		}

}


}

///Decode a normal pdb frame 
void decode_pdb_frame(struct expbuf *buf,struct system *in,int justgetinfo)
{
int quit=FALSE;
int i=0;
int ii;
char in1[20];
int last_mol=1;
int atoms=0;
int mol=1;	// internal pdb file number.
char buffer[1000];
int serial;
char space1;
char name[10];
char var;
char residname[10];
char space2;
char chainID;
char space3;
char space4;
char space5;
char space7;
char xstr[20];
char ystr[20];
char zstr[20];


system_free(in);


//printf("length now 2 = -%d-\n",in->length);
//getchar();
if (in->length<=0)
{
	//printf("\n\nalloc\n\n");
	mol_expand_memory(in,0);
}

ii=-1;
atoms=0;
do
{


		//fgets(buffer, 1000, file);
		//printf("aa\n");
		expbuf_gets(buffer,buf);
		//printf("a: %s-a-\n",buffer);
		//printf("b: %s-b-\n",buf->buf);
		//getchar();
		
		if ((mycmp4(buffer,"ATOM"))||(mycmp(buffer,"HETATM")))
		{
			//printf("atom used to process\n");

			if (justgetinfo==TRUE)
			{
				quit=TRUE;
				break;
			}
			//This is not ideal but pre read it to figire out what is going on
			sscanf(buffer,"%6s%5d%c%4s%c%3s%c%c%4d%c%c%c%c%8c%8c%8c",in1,&serial,&space1,name,&var,residname,&space2,&chainID,&mol,&space7,&space3,&space4,&space5,xstr,ystr,zstr);

			if ((atoms==0)&&(last_mol==1))  last_mol=mol;


			if (last_mol!=mol)
			{
				//printf("%d %d\n",in->length,i);
				if (in->length<=(i+1))
				{
					//printf("alloc\n");
					mol_expand_memory(in,1);
				}else
				{
					//printf("No alloc\n");
				}
				ii=0;
				i++;
			}else
			{
				//printf("%d\n",in->mols[i].atoms);
				if (in->mols[i].atoms<=(ii+1))
				{
					//printf("alloc2\n");
					mol_resize(&(in->mols[i]),1);
				}else
				{
					//printf("no alloc2\n");
				}
				ii++;
			}

			last_mol=mol;
			atoms++;

			set_null(in->mols[i].element[ii],4);
			sscanf(buffer,"%6s%5d%c%2c%c%c%c%3s%c%c%4d%c%c%c%c%8c%8c%8c",in1,&serial,&space1,in->mols[i].element[ii],&(in->mols[i].dist[ii]),&(in->mols[i].branch[ii]),&var,in->mols[i].resid,&space2,&(in->mols[i].chainid[ii]),&mol,&space7,&space3,&space4,&space5,xstr,ystr,zstr);

			xstr[8]=0;
			ystr[8]=0;
			zstr[8]=0;

			sscanf(xstr,"%lf",&(in->mols[i].atom[ii].x));
			sscanf(ystr,"%lf",&(in->mols[i].atom[ii].y));
			sscanf(zstr,"%lf",&(in->mols[i].atom[ii].z));


			remove_white_space(in->mols[i].element[ii]);

			in->mols[i].vel[ii].x=-100.0;
			in->mols[i].vel[ii].y=-100.0;
			in->mols[i].vel[ii].z=-100.0;


		}else
		if (mycmp(buffer,"CRYST1"))
		{
			//printf("crystal used to process\n");
			sscanf(buffer,"%s %lf %lf %lf",in1,&(in->size.x),&(in->size.y),&(in->size.z));
		}else
		if (mycmp(buffer,"TITLE"))
		{
			//TITLE     Protein t=   0.00000
			int pos=util_find_end_of(buffer,"t=");
			if (pos!=-1)
			{
				sscanf((buffer+pos+2),"%lf",&(in->time));
			}
		}else
		if (mycmp(buffer,"MODEL"))
		{
			//printf("Model used\n");
			//TITLE     Protein t=   0.00000
			int pos=util_find_end_of(buffer,"MODEL");
			if (pos!=-1)
			{
				sscanf((buffer+strlen("MODEL")),"%d",&(in->model));

			}
			//printf("%s -%d-",buffer,pos);
			//printf("%d\n",in->model);
		}
		else
		if (mycmp(buffer,"END"))
		{
			quit=TRUE;
		}else
		{
			//printf("Not processed\n");
		}

		if (expbuf_eof(buf)==TRUE)
		{
			quit=TRUE;
		}

		//printf("atoms= %d\n",in->mols[0].atoms);
		//getchar();
}while (quit==FALSE);

}

///Load a single pdb frame from a file and then close the file.
void load_pdb(char * file_name,struct system *in)
{
FILE *file;
if((file=fopen(file_name, "rb")) == NULL)
{
	printf("Cannot open file for read.\n");
	exit(1);
}

load_pdb_frame(file,in,FALSE);

fclose(file);
}

///Load a compressed pdb frame from a file then close it.
void load_pdb_comp(char * file_name,struct system *in)
{
FILE *file;
if((file=fopen(file_name, "rb")) == NULL)
{
	printf("Cannot open file for read.\n");
	exit(1);
}

load_pdb_frame_comp(file,in,FALSE);

fclose(file);
}

int isnumber(char in)
{
if (in=='0') return TRUE;
if (in=='1') return TRUE;
if (in=='2') return TRUE;
if (in=='3') return TRUE;
if (in=='4') return TRUE;
if (in=='5') return TRUE;
if (in=='6') return TRUE;
if (in=='7') return TRUE;
if (in=='8') return TRUE;
if (in=='9') return TRUE;
return FALSE;
}

///Load a gro (gromacs) file.  This probably needs rewriting.  
void load_gro(char * name,struct system *in)
{
char buf[500];
char temp[500];
int el=0;
char b2[20];
char b3[20];
double last_el=0;
int atom;
double p1;
double p2;
double p3;
int atoms=0;
double v1;
double v2;
double v3;
int total_atoms=0;

FILE *file;
if((file=fopen(name, "r")) == NULL)
{
	printf("Cannot open file for read.\n");
	exit(1);
}

fgets(buf, 500, file);
fgets(buf, 500, file);
fgets(buf, 500, file);
int delta=0;
int i;
for (i=1;i<strlen(buf);i++)
{
if ((buf[i-1]!=' ')&&(buf[i]==' ')) delta++;
}
//printf("delta= %d\n",delta);
//getchar();
rewind(file);

fgets(buf, 500, file);
int molcount=0;
fscanf(file, "%d",&total_atoms);
//printf("length = %d\n",total_atoms);
//getchar();
int count=0;
do
{
	if (delta>5)
	{
		fscanf(file,"%5d%5c%5c%5d%lf%lf%lf%lf%lf%lf",&el,b2,b3,&atom,&p1,&p2,&p3,&v1,&v2,&v3);
	}else
	{
		fscanf(file,"%5d%5c%5c%5d%lf%lf%lf",&el,b2,b3,&atom,&p1,&p2,&p3);
		v1=-100.0;
		v2=-100.0;
		v3=-100.0;
	}
	b2[5]=0;
	b3[5]=0;
	remove_white_space(b2);
	remove_white_space(b3);
/*if (count>9996)
{
	printf("%d|%s|%s|%d|%f|%lf|%8.3lf %8.4lf %8.4lf %8.4lf\n",el,b2,b3,atom,p1,p2,p3,v1,v2,v3);
	getchar();
}*/
	if (count!=0)
	if (el!=last_el)
	{
		mol_expand_memory(in,atoms);
		atoms=0;

	}



	last_el=el;
	count++;

	atoms++;

}while (count<total_atoms);

mol_expand_memory(in,atoms);

fscanf(file,"%lf %lf %lf",&in->size.x,&in->size.y,&in->size.z);

in->size.x*=10.0;
in->size.y*=10.0;
in->size.z*=10.0;

rewind(file);

fgets(buf, 500, file);
fscanf(file, "%d",&total_atoms);
//printf("length = %d\n",total_atoms);

molcount=0;
count=0;
atoms=0;
do
{
	set_null(b2,20);
	set_null(b3,20);
	if (delta>5)
	{
		fscanf(file,"%5d%5c%5c%5d%lf%lf%lf%lf%lf%lf",&el,b2,b3,&atom,&p1,&p2,&p3,&v1,&v2,&v3);
	}else
	{
		fscanf(file,"%5d%5c%5c%5d%lf%lf%lf",&el,b2,b3,&atom,&p1,&p2,&p3);
		v1=-100.0;
		v2=-100.0;
		v3=-100.0;
	}

	remove_white_space(b2);
	remove_white_space(b3);

//printf("here1!\n");
	if (count!=0)
	if (el!=last_el)
	{
		//mol_expand_memory(mols,nmols,atoms);
		//strcpy(last,b1);
		atoms=0;
		molcount++;
	}
//printf("here2 %d %d!\n",molcount,atoms);

//getchar();


	in->mols[molcount].atom[atoms].x=p1*10.0;
	in->mols[molcount].atom[atoms].y=p2*10.0;
	in->mols[molcount].atom[atoms].z=p3*10.0;
//printf("here3!\n");

	in->mols[molcount].vel[atoms].x=v1;
	in->mols[molcount].vel[atoms].y=v2;
	in->mols[molcount].vel[atoms].z=v3;
	//printf("one '%s'\n",b3);
	int j=0;
	int max=3-strlen(b3);

	if (strlen(b3)==1)
	{
	strcpy(temp," ");
	strcat(temp,b3);
	strcpy(b3,temp);
	}
//I just put this in?
	if (strlen(b3)==2)
	{
	strcpy(temp," ");
	strcat(temp,b3);
	strcpy(b3,temp);
	}


	if ((strlen(b3)==3)&&(isnumber(b3[strlen(b3)-1])))
	{
	strcpy(temp," ");
	strcat(temp,b3);
	strcpy(b3,temp);
	}


	if (max<0)
	{
		printf("%s Error!\n",b3);
		exit(0);
	}

	for (j=0;j<max;j++)
	{
		strcat(b3," ");
	}

	in->mols[molcount].element[atoms][0]=b3[0];
	in->mols[molcount].element[atoms][1]=b3[1];
	in->mols[molcount].element[atoms][2]=0;
	remove_white_space(in->mols[molcount].element[atoms]);
	//printf("'%s'\n",b3);
	//getchar();
	in->mols[molcount].dist[atoms]=b3[2];
	in->mols[molcount].branch[atoms]=b3[3];

	in->mols[molcount].chainid[atoms]=' ';

	strcpy(in->mols[molcount].resid,b2);
//printf("here4!\n");

	last_el=el;
	count++;

	atoms++;

}while (count<total_atoms);



fclose(file);
}

///Load an xyz file, like the ones from gaussian.
void load_xyz(struct system *in,char *input)
{
}
