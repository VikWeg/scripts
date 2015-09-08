#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <functional>
#include <random>
#include <ctime>
#include <direct.h>
#include <string>
#include <cmath>
#include <gnuplot_i.hpp>
#include "nifti1.h"

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352
nifti_1_header hdr;

float* data;
float* L1data;
float* L2data;
float* L3data;
float* mask;

std::ofstream par_file;
std::ofstream Efile;
FILE* fiber_file;

char buffer[80];

struct vertex
{
	float x;
	float y;
	float z;

	float pos_x;
	float pos_y;
	float pos_z;

	float T0;
	float T1;
	float T2;
	float T3;
	float T4;
	float T5;

	int* c;

	float Emin;
	float Emax;
	float delta_E;

	int cc;
	int nn;
	int sig;

	vertex** n;
};

vertex**** ensemble;
int*** snum;

int*** wmask;
int*** surf_mask;

long w_vox_num;
long surf_vox_num;
long scount;

const float eps = 10e-7;

// **** PARAMETERS ********************************************************/
/**/	std::string version("Linked List CPU");	/**/
/**/	std::string subject ("1159T");									/**/
/**/	int cube_size[3] = { 9, 9,9 };									/**/
/**/	int vox_origin[4] = { 0, 48, 64, 64 };							/**/
/**/																	/**/
/**/	float cutoff = 1.5;												/**/
/**/																	/**/
/**/	float Ti = 0.05;												/**/
/**/	float Tf = 0.00005;												/**/
/**/	float etha = 0.9;												/**/
/**/	long tsteps_tot = ceilf(logf(Tf / Ti) / logf(etha));			/**/
/**/																	/**/
/**/	int nx = 30;													/**/
/**/	int S = 1;														/**/
/**/	float delta_x = 0.1;											/**/
/**/																	/**/
/**/	char* wc_str = "1/T";												/**/
/**/	char* wx_str = "1";												/**/
/**/																	/**/
/**/	float wc(float T) { return 1/T; }								/**/
/**/	float wx(float T) { return 1; }	
/**/
/**/	float wint(float cos) { return (1 + cos) / (1-cos + 0.0000001); }
/**/	char* wint_str = "(1+cos)/2";

		char* wdata_str = "(E1 + E2)/2";
/**/
/**/	std::string comment("");
/*************************************************************************/

// **** RANDOM ****

std::default_random_engine generate;
std::uniform_int_distribution<int> u_i(0, cube_size[0] - 1);
std::uniform_int_distribution<int> u_j(0, cube_size[1] - 1);
std::uniform_int_distribution<int> u_k(0, cube_size[2] - 1);

// **** TIMING ****

float CPU_freq = 4*pow(10., 9.);

typedef union
{
	__int64 int64;
	struct { __int32 lo, hi; } int32;
} tsc_counter;

#define RDTSC(cpu_c) \
{ __asm rdtsc \
	__asm mov(cpu_c).int32.lo, eax \
	__asm mov(cpu_c).int32.hi, edx \
}

tsc_counter start, stop, start_all, stop_all;