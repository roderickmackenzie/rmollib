/** @file rnd.h
	@brief Header file for rnd.c
*/
#ifndef rndh
#define rndh
//Random number generator
int rnd_init ();
void rnd_free();
double get_random_number(double in);
void id_rand(unsigned char *buffer,int length);
#endif
