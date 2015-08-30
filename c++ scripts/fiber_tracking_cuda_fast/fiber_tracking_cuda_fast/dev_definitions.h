#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cuda.h"
#include "cuda_runtime_api.h"
#include "device_functions.h"

#include <curand.h>
#include <curand_kernel.h>

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
int dim[4] = { 0, 0, 0, 0 };

__device__ float* dev_pixdim;

float* data;
float* L1data;
float* L2data;
float* L3data;
float* mask;

std::ofstream par_file;
std::ofstream Efile;
FILE* fiber_file;

char buffer[80];

// ========= Vertex Arrays =============

//read & write
float* x;
float* y;
float* z;

int* cc;
int* c; //accessed via n_id;

// read-only
float* pos_x;
float* pos_y;
float* pos_z;

float* T0;
float* T1;
float* T2;
float* T3;
float* T4;
float* T5;

float* Emin;
float* Emax;
float* delta_E;

int* nc;
int* sig;

int* n_id;
int* n; //accessed via n_id;

// ======================

//*** read & write ***

// dev_in
float* dev_in_x;
float* dev_in_y;
float* dev_in_z;

int* dev_in_cc;
int* dev_in_c;

//dev_out
float* dev_out_x;
float* dev_out_y;
float* dev_out_z;

int* dev_out_cc;
int* dev_out_c;

//*** read-only ***
float* dev_pos_x;
float* dev_pos_y;
float* dev_pos_z;

float* dev_T0;
float* dev_T1;
float* dev_T2;
float* dev_T3;
float* dev_T4;
float* dev_T5;

float* dev_Emin;
float* dev_Emax;
float* dev_delta_E;

int* dev_nc;
int* dev_sig;

int* dev_n_id;
int* dev_n;

// ======================
int* snum;

int* wmask;
int* surf_mask;

long w_vox_num;
long surf_vox_num;
long scount;

long n_num;

struct parameters
{
	int nx;
	float delta_x;
	int scount;
	float pix_dim_x;
	float pix_dim_y;
	float pix_dim_z;
};

parameters* host_params;
parameters* dev_params;

// How does voxel size affect cost function?
// What if all tensors are isotropic -> no solution?
// Make sure temperature scaling is calibrated, so that energy does not diverge at T -> 0.
// use multiple of 32 for number of threads
// use struct of arrays (coalesced memory)
// keep parameters (nx, delta_x, T, scount...) in const memory!
// try to avoid unnecessary mc_x evaluations of E0 after rejection.
// be aware that there is a conflict when updating c -> can be updated by different threads (leaving it for now).


// **** PARAMETERS ************************************************/
/**/															/**/
/**/	std::string subject("1159T");							/**/
/**/	int cube_size[3] = { 9, 9, 9 };							/**/
/**/	int vox_origin[4] = { 0, 48, 64, 64 };					/**/
/**/															/**/
/**/	float cutoff = 1.5;										/**/
/**/															/**/
/**/	float Ti = 0.025;										/**/
/**/	float Tf = 0.0005;										/**/
/**/	float etha = 0.5;										/**/
/**/	int tsteps_tot = (int)ceilf(logf(Tf / Ti) / logf(etha));	/**/
/**/															/**/
/**/	int nx = 30;								/**/
		__device__ int* dev_nx;
/**/	int sweeps = 2;												/**/
/**/	float delta_x = 0.1;									/**/
		__device__ float* dev_delta_x;
/**/															/**/
/**/	char* wc_str = "1 / T";									/**/
/**/	char* wx_str = "1";						/**/
/**/															/**/
/**/	float wc(float T) { return 1 / T; }					/**/
/**/	float wx(float T) { return 1; }
/**/	__device__ float dev_wc(float T) { return 1 / T; }					/**/
/**/	__device__ float dev_wx(float T) { return 1; }

/**/
/**/	float wint(float cos) { return 0.5*(1 + cos); }
/**/	__device__ float dev_wint(float cos) { return 0.5*(1 + cos); }
/**/	char* wint_str = "(1+cos)/2";

char* wdata_str = "(E - Emin)/(Emax-Emin)";
/**/
/**/	std::string comment("");
/******************************************************************/

// **** TIMING ****

float CPU_freq = 4 * pow(10., 9.);

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