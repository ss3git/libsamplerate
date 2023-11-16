/*
** Copyright (c) 2002-2016, Erik de Castro Lopo <erikd@mega-nerd.com>
** All rights reserved.
**
** This code is released under 2-clause BSD license. Please see the
** file at : https://github.com/libsndfile/libsamplerate/blob/master/COPYING
*/

#define	ABS(a)			(((a) < 0) ? - (a) : (a))

#ifndef MAX
#define	MAX(a,b)		(((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define	MIN(a,b)		(((a) < (b)) ? (a) : (b))
#endif

#define	ARRAY_LEN(x)	((int) (sizeof (x) / sizeof ((x) [0])))

void gen_windowed_sines (int freq_count, const double *freqs, double max, float *output, int output_len) ;
void gen_windowed_sines_S32 (int freq_count, const double *freqs, double max, int *output, int output_len) ;
void gen_windowed_sines_D64 (int freq_count, const double *freqs, double max, double *output, int output_len) ;

void save_oct_float (char *filename, float *input, int in_len, float *output, int out_len) ;
void save_oct_double (char *filename, double *input, int in_len, double *output, int out_len) ;

void interleave_data (const float *in, float *out, int frames, int channels) ;

void deinterleave_data (const float *in, float *out, int frames, int channels) ;

void reverse_data (float *data, int datalen) ;

double calculate_snr (float *data, int len, int expected_peaks) ;
double calculate_snr_S32 (int *data, int len, int expected_peaks) ;
double calculate_snr_D64 (double *data, int len, int expected_peaks) ;

const char * get_cpu_name (void) ;

#define ASSERT(condition) \
	if (!(condition)) \
	{	printf ("Condition failed on Line %d : %s\n\n", __LINE__, #condition) ; \
		exit (1) ; \
	    } ;


enum { // these are only for test suites.
	SRC_S32_START = 32,
#ifdef SUPPORT_S32_INTERFACE
	SRC_SINC_BEST_QUALITY_S32 = SRC_S32_START,
#endif
	SRC_D64_START = 64,
#ifdef SUPPORT_S32_INTERFACE
	SRC_SINC_BEST_QUALITY_D64 = SRC_D64_START,
#endif
};
