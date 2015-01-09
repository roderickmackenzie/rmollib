/** @file save.h
	@brief Header file for save.c
*/
#ifndef saveh
#define saveh
void save_pdb_frame_comp(FILE *out,double shift,struct system *in);
void build_pdb_frame(struct expbuf *out,double shift,struct system *in);
void save_file(char * name,double shift,struct system *mysystem);
void dump_com_file(struct mol *my_vec1,struct mol *my_vec2, char *name,struct vec *cog,char *gaus_base_file);
void dump_box_as_pdb(char *name,double shift,struct system *in);
void dump_pair(struct mol *my_vec1,struct mol *my_vec2,char *name,struct vec *shift);
void dump_box_as_gro(struct mol *myarray,int length, char *name,struct system *mysystem);
void save_pdb_frame(FILE *out,double shift,struct system *in);
void save_pdb_comp(char *name,double shift,struct system *in);
void build_z_frame(unsigned char *out,double shift,struct system *in);
#endif
