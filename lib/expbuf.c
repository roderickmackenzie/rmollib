/** @file expbuf.c
	@brief A set of functions to accept any amount of string data from a buffer of any size and then read it out line by line.  It is very useful when reading a file containing strings, storing it in memory and then reading it back line by line later.  
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include "true.h" 
#include "expbuf.h"
/**Initialize the expbuf structure.  This function must be called before using expbuf 
@param in A pointer to the buffer to be initialized. 
*/
void expbuf_init(struct expbuf * in)
{
	in->len=0;
	in->pos=0;
	in->max=10;
	in->buf=NULL;
	in->buf=(char *)realloc(in->buf,sizeof(char)*in->max);
	strcpy(in->buf,"");
}

/**Free all memory associated with the expbuf structure.  This function must be called before using expbuf 
@param in A pointer to the buffer to be freed. 
*/
void expbuf_free(struct expbuf * in)
{
	in->len=0;
	in->max=0;
	in->pos=0;
	free(in->buf);
}


/**Add a string to the expbuf buffer
@param in A pointer to the buffer to add data to.
@param add The data to add to the buffer 
*/
void expbuf_add(struct expbuf * in,char *add)
{
int len=strlen(add);

if ((in->len+len)>=in->max)
{
	in->max+=10240;
	in->buf=(char *)realloc(in->buf,sizeof(char)*in->max);
//printf("grown\n");
}else
{
//printf("not grown\n");

}

strcpy((in->buf+in->len),add);
in->len+=len;
}


/**Check if the end of the buffer has been reached.
@param in A pointer to the buffer to check.
*/
int expbuf_eof(struct expbuf * in)
{
if (in->pos>=in->len)
{
	return TRUE;
}
return FALSE;
}

/**Get a string from the buffer
@param in A pointer to the buffer to extract data from.
@param out a pointer to the buffer to which the string should be transfered.
*/
void expbuf_gets(char *out,struct expbuf * in)
{
char * outbuild=out;
char * inbuild=(in->buf+in->pos);
int pos=in->pos;
do
{
*(outbuild++)=*(inbuild++);
pos++;
	if ((*inbuild)=='\n')
	{
		inbuild++;
		pos++;
		break;
	}

 if (pos>=in->len) break;

}while((*inbuild)!=0);
*(outbuild++)=0;

in->pos=pos;
}
