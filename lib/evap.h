/** @file evap.h
	@brief Header file for evap.c
*/
#ifndef evaph
#define evaph
//evap
void evap_edge(struct system * in,char *name,double percent);
void evap_min(struct system * in,char *name,double percent);
void evap(struct system * in,char *name,double percent);

#endif
