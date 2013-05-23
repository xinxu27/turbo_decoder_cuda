#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

typedef unsigned char BYTE;
typedef int INT;
typedef unsigned int UINT;
typedef int BOOL;

#define INIFINITY  1E+10

void TurboCodingInit();
void TurboEnCoding(int *source, int *coded_source, int source_length);
void module(int * a,double * outi,double * outq,int N, int modu_index);
void AWGN(double *send, double *r, double sigma, int totallength);
void demodule(double *symbol_i, double *symbol_q, int symbol_len,float* out,double Kf,int modu_index);

void TurboDecoding(float *flow_for_decode, int *flow_decoded,int flow_length);

void TurboCodingRelease();

void countErrors(BYTE *m, BYTE * mhat, UINT * bitsError, UINT * frameError, UINT iter);
__global__ void demultiplex(float * stream, float * msg, float * parity);
