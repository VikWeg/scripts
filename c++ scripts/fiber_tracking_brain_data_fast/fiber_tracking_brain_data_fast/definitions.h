#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <functional>
#include <random>
#include "nifti1.h"

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352
nifti_1_header hdr;
float* data;
float* L1data;
float* L2data;
float* L3data;
float* mask;

std::ofstream fiber_file;

std::default_random_engine generate;

struct vertex
{
	float x;
	float y;
	float z;

	int pos_x;
	int pos_y;
	int pos_z;

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

int size[3] = { 9, 9, 9 };
int center[4] = { 0, 48, 64, 64 };

int*** wmask;
int*** surf_mask;

long w_vox_num;
long surf_vox_num;

std::uniform_int_distribution<int> u_i(0, size[0] - 1);
std::uniform_int_distribution<int> u_j(0, size[1] - 1);
std::uniform_int_distribution<int> u_k(0, size[2] - 1);

float cutoff = 1.5;

float Ti = 0.3;
float Tf = 0.002;
float etha = 0.9;

int nx = 30;
//int nc = 50;
int S = 5;

float delta_x = 0.1;
float wc(float T) {return 1./T;}
float wx(float T) {return 1;}


// **** TIMING ****

float CPU_freq = 4.*pow(10., 9.);

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

tsc_counter start, stop;