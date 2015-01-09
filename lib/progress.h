/** @file progress.h
	@brief Header file for progress.c
*/
#ifndef progress_h
#define progress_h
#define fg_reset	0
#define fg_red		31
#define fg_green	32
#define fg_yellow	33
#define fg_blue		34
#define fg_purple	35

void progress_clear();
void set_progress_multi_line_text(char *in);
void set_progress_multi_line();
void text_progress(double percent);
void text_progress_finish();
void text_progress_start(char *in);
void set_porgress_max(int in);
void set_porgress_noreset();
void set_progress_colored();
void set_porgress_nospin();
void set_porgress_color(int in);

#endif
