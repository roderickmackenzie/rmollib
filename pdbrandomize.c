/** @file pdbrandomize.c
	@brief This program reads in a pdb or gro file and randomly changes the order of the molecules with in the file.  This is to make it look nicer when pymol plots the images.  It destroys the rainbow effect so the individual molecules stand out.
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


int main(int argc, char* argv[])
{
//mtrace();
struct system mainsystem;
system_init(&mainsystem);

if (scanarg( argv,argc,"--in")==FALSE)
{
	printf("Need in\n");
}

if (scanarg( argv,argc,"--out")==FALSE)
{
	printf("need --out!\n");
	exit(0);
}

load_file(get_arg_plusone( argv,argc,"--in"),&mainsystem);
system_post_load_fix_up(&mainsystem);
int i=0;
int rnd;
unsigned char temp[sizeof(struct mol)];
for (i=0;i<mainsystem.length;i++)
{
	rnd=(int)get_random_number((double)(mainsystem.length-1));
	memcpy(&temp,&(mainsystem.mols[i]),sizeof(struct mol));
	memcpy(&mainsystem.mols[i],&(mainsystem.mols[rnd]),sizeof(struct mol));
	memcpy(&(mainsystem.mols[rnd]),&temp,sizeof(struct mol));
}

save_file(get_arg_plusone( argv,argc,"--out"),0.0,&mainsystem);


system_free(&mainsystem);

return 0;

}
