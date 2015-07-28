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
	float* x;
	int* c;
	int cc;
	float* T;
	int* pos;
	int sig;
	float Emin;
	float Emax;
	int nn;
	vertex** n;
};

vertex**** ensemble;
int*** snum = new int**;

int size[3] = { 20, 20, 20 };
int center[4] = { 0, 93, 59, 29 };

int*** wmask;
int*** surf_mask;

long w_vox_num;
long surf_vox_num;

std::uniform_int_distribution<int> u_i(0, size[0] - 1);
std::uniform_int_distribution<int> u_j(0, size[1] - 1);
std::uniform_int_distribution<int> u_k(0, size[2] - 1);

float cutoff = 1.5;
float cosmax = 0;

float Ti = 1.5;
float Tf = 0.001;
float etha = 0.8;

int nx = 30;
int nc = 30;
int S = 5;

float delta_x = 0.1;
float wc(float T) {return 1./T;}
float wx(float T) {return 1;}

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