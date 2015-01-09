/** @file random.c
	@brief A random number generator.
*/
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
//Rods random number generator
#include <stdio.h>
#include <gsl/gsl_rng.h>

const gsl_rng_type *gsl_rng_default;
unsigned long int gsl_rng_default_seed;

static gsl_rng * r;

///Initialize the random number generator
int rnd_init ()
{
const gsl_rng_type * T;


gsl_rng_env_setup();

T = gsl_rng_taus;//gsl_rng_ranlux389;

//gsl_rng_default_seed=time(NULL);
r = gsl_rng_alloc (T);
//printf("here------------------------------\n");
FILE *infile;
infile =fopen("/dev/urandom", "rb");
long int start;
fread(&start, sizeof(long int), 1, infile);
fclose(infile);

gsl_rng_set (r, start);


return 0;
}

///Free the random number generator
void rnd_free()
{
gsl_rng_free (r);
}

///generate a random number
double get_random_number(double in)
{
	return gsl_rng_uniform (r)*in;
}

///generate a buffer full of random chars
void id_rand(unsigned char *buffer,int length)
{
int i;
double u;
	for (i = 0; i < length; i++) 
	{
		u = gsl_rng_uniform (r);
		buffer[i]=(unsigned char)(u*256);
	}
}
