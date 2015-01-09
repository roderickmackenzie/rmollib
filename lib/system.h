/** @file system.h
	@brief Header file for system.c
*/
/**The system header files provides access to molecular systems comprised of molecules
*/
#ifndef systemh
#define systemh

#define _FILE_OFFSET_BITS 64



///Structure defining a box minima and maxima
struct box
{
double xmin;
double xmax;
double ymin;
double ymax;
double zmin;
double zmax;
};


///Structure to hold an entire molecular system
struct system
{
struct mol *mols;
int length;
struct vec size;
struct box systembox;
struct vec fsize;
struct box funcbox;
struct vec cog;
double time;
double r;
int model;
};
void system_insert_mol(struct system *out,struct mol *in,struct vec *shift);
int system_count_resid(struct system * in,char *name);
void system_set_box_now(struct system *in);
void system_cal_cog_n(struct system *in,int n);
void system_inv_del_mols(struct system *in , char *resid );
void system_del_mols(struct system *in,char *resid);
void system_find_resids(struct system * in);
void system_post_load_fix_up(struct system *in);
double system_diff(struct system *one,struct system *two,char *resid);
void system_get_min(struct system *in,struct vec *ret);
void system_init(struct system* in);
void system_set_group(struct system *mysystem,int len);
void rot_system(struct system *mysystem,struct vec *unit,struct vec *base, double ang);
double system_cal_density(struct system *in);
void system_cal_step(struct vec *step,struct system *in);
void add_to_system(struct system *out,struct system *in,struct vec *shift);
void system_shift(struct system *in,struct vec *delta);
void system_rotx(struct system *in,double ang);
void system_roty(struct system *in,double ang);
void system_rotz(struct system *in,double ang);
void system_free(struct system *in);
void system_set_box(struct system *in);
void system_cal_cog(struct system *in);
void system_rescale(struct system *in);
void system_truncate_v(struct system *in,double v);
void system_expand(struct system *out,struct system *in,double size_x,double size_y,double size_z, int tidy);
struct mol * mol_expand_memory(struct system * in,int atoms);
void system_cpy(struct system *out,struct system *in);
void system_cal_r(struct system *in);
#endif
