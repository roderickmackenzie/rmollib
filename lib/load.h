/** @file load.h
	@brief Header file for load.c
*/
#ifndef load_h
#define load_h
//Loading and saving routines
void load_file(char * name,struct system *in);
void load_pdb(char * file_name,struct system *in);
void load_gro(char * file_name,struct system *in);
void load_xyz(struct system *in,char *input);
void decode_pdb_frame(struct expbuf *buf,struct system *in,int justgetinfo);
int load_pdb_frame(FILE *file,struct system *in,int justgetinfo);
void load_pdb_frame_comp(FILE *file,struct system *in,int justgetinfo);
void load_pdb_comp(char * file_name,struct system *in);
void decode_z_frame(unsigned char *buf,struct system *in,double time);

#endif

