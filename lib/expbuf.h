/** @file expbuf.h
	@brief Header file for expbuf.c
*/

#ifndef expbufh
#define expbufh
struct expbuf
{
///The lenght of the data stored in the buffer
long int len;
///The maximum length of the buffer, this is bigger than len.  This is because allocation happens in blocks for maximum speed.
long int max;
///Current buffer read position
long int pos;
///Pointer to the buffer
char *buf;
};

void expbuf_init(struct expbuf * in);
void expbuf_free(struct expbuf * in);
void expbuf_add(struct expbuf * in,char *add);
int expbuf_eof(struct expbuf * in);
void expbuf_gets(char *out,struct expbuf * in);

#endif
