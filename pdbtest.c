/** @file pdbtest.c
	@brief Program to the function of the pdb and gro read/write.
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
	if (get_arg_plusone_pos( argv,argc,"--out")!=-1)
	{
	}else
	{
		printf("need --out!\n");
		exit(0);
	}

	if (get_arg_plusone_pos( argv,argc,"--in")!=-1)
	{
	}else
	{
		printf("need --in!\n");
		exit(0);
	}

	load_file(get_arg_plusone( argv,argc,"--in"),&mainsystem);
	//printf("length now = %d\n",mainsystem.length);
	//load_file(get_arg_plusone( argv,argc,"--in"),&mainsystem);
	//string_endcap(&mainsystem,TRUE,TRUE);
	save_file(get_arg_plusone( argv,argc,"--out"),0.0,&mainsystem);
	system_free(&mainsystem);


system_free(&mainsystem);

return 0;

}
