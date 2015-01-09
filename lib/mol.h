/** @file mol.h
	@brief Header file for mol.c
*/
#ifndef molh
#define molh
#include <stdint.h>
#define		pi	3.14159265358979323846
#include "true.h" 
#include "vec.h"
#include "expbuf.h"
#include "rnd.h"
#include "molops.h"
#include "system.h"
#include "util.h"
#include "evap.h"
#include "load.h"
#include "save.h"
#define _FILE_OFFSET_BITS 64

///Header for pdbz file
struct comp_header
{
///Compressed length of the frame
uint64_t len;
///Original length of the file
uint64_t origlen;
///Time of the frame	
double time;
///x size of simulation box
double xsize;
///y size of simulation box
double ysize;
///z size of simulation box
double zsize;
};


double ret_ram(char *in);
double get_max(double *in,int length);

#endif
