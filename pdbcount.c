/** @file pdbcount.c
	@brief Program to count the number of residues in a file.
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
printf("need --in\n");
exit(0);
}


load_file(get_arg_plusone( argv,argc,"--in"),&mainsystem);
system_post_load_fix_up(&mainsystem);


save_file("./save.pdb",0.0,&mainsystem);
if (get_arg_plusone_pos( argv,argc,"--resid")!=-1)
{
}else
{
	printf("Tell me what residule name you want!\n");
	exit(0);
}

int number=system_count_resid(&mainsystem,get_arg_plusone( argv,argc,"--resid"));
printf("%d\n",number);

system_free(&mainsystem);

return 0;

}
