/** @file util.h
	@brief Header file for util.c
*/
#ifndef utilh
#define utilh
//Utility routines
void copyright();
void platofrm_test();
int scanarg( char* in[],int count,char * find);
int get_arg_plusone_pos( char* in[],int count,char * find);
char * get_arg_plusone( char* in[],int count,char * find);
void right_align(char *out, int space,char *in);
void remove_white_space(char *in);
double ret_ram(char *in);
int mycmp4(char *buf,char *in);
int mycmp(char *buf,char *in);
#endif
