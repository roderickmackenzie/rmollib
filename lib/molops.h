/** @file molops.h
	@brief Header file for molops.c
*/
#ifndef molopsh
#define molopsh

///Structure to hold molecule
struct mol
{
double theta;		///delete?
double rot;		///delete?
double alpha;		///delete?
double cos_;		///delete?
int group;		///group number
struct vec n;
int *chums;		///list of neighbors
int nchums;
struct vec cog;		///Center of gravity
struct vec cog_n;	///Center of gravity for first n points
struct vec * atom;
struct vec * vel;	///Velocity
char *chainid;
char *dist;
char *branch;
char **element;
int atoms;
char *resid;
double r;
double mass;
int del;
};

//Manipulation on the molecular level
int mol_check_bond(struct mol *my_mol1,int chum1,struct mol *my_mol2,int chum2);
double mol_diff(struct mol *one,struct mol *two);
void mol_rot(struct mol *mymol,struct vec *unit,struct vec *base, double ang);
void mol_zero(struct mol *my_mol,int n);
void mol_rotx(struct mol *in,double ang);
void mol_roty(struct mol *in,double ang);
void mol_rotz(struct mol *in,double ang);
void mol_init(struct mol *my_mol,int atoms);
void mol_cal_cog(struct mol *my_mol);
void mol_shif(struct mol *in,struct vec *delta);
void mol_cpy(struct mol *my_mol1,struct mol *my_mol2);
void mol_free(struct mol *my_mol);
void mol_add_atom_to(struct mol *my_mol,char * name,struct vec *pos);
int mol_cal_spin(struct mol *in1,struct mol *in2);
void mol_cal_r(struct mol *my_mol);
void mol_cal_step(struct mol *my_mol);
void mol_get_ang(struct mol *my_mol,struct mol *my_mol2);
void mol_set_mass(struct mol *my_mol);
void mol_resize(struct mol *my_mol,int atoms);
void mol_cal_r(struct mol *my_mol);
void mol_get_ang(struct mol *my_mol,struct mol *my_mol2);
void mol_set_mass(struct mol *my_mol);
void mol_shift(struct mol *in,struct vec *delta);
void mol_cal_cog_n(struct mol *my_mol,int n);
#endif
