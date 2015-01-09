/** @file pdbevap.c
	@brief This program reads a single frame of a pdb file and removes at random a proportion of a desired molecule.  It is designed for solvent removal i.e. evaporation.  It has three modes of operation, removal of solvent at random.  Removal of solvent from the edges of a simulation cell.  Removal of solvent that is furthest away from its neighbors.    
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
rnd_init ();
struct system mainsystem;
system_init(&mainsystem);

double per=0.0;

if (scanarg( argv,argc,"--in")==FALSE)
{
	printf("Need --in\n");
	exit(0);
}

	if (get_arg_plusone_pos( argv,argc,"--remove")!=-1)
	{
		per=atof(get_arg_plusone( argv,argc,"--remove"));
	}else
	{
		printf("Tell me how many molecules to remove!\n");
		exit(0);
	}

	if (get_arg_plusone_pos( argv,argc,"--out")!=-1)
	{
	}else
	{
		printf("Tell me where you want me to dump it!\n");
		exit(0);
	}

	if (get_arg_plusone_pos( argv,argc,"--sol")!=-1)
	{
	}else
	{
		printf("Tell me which one is the solvent --sol!\n");
		exit(0);
	}


	load_file(get_arg_plusone( argv,argc,"--in"),&mainsystem);
	system_post_load_fix_up(&mainsystem);

	if (scanarg( argv,argc,"--takemin"))
	{
		evap_min(&mainsystem,get_arg_plusone( argv,argc,"--sol"),per);
	}else
	if (scanarg( argv,argc,"--takeedge"))
	{
		printf("doing edge\n");
		evap_edge(&mainsystem,get_arg_plusone( argv,argc,"--sol"),per);
	}else
	{
		evap(&mainsystem,get_arg_plusone( argv,argc,"--sol"),per);
	}
	

	save_file(get_arg_plusone( argv,argc,"--out"),0.0,&mainsystem);


system_free(&mainsystem);
rnd_free ();
return 0;

}
