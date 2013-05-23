#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

void TurboCodingInit();
void TurboEnCoding(int *source, int *coded_source, int source_length);
void module(int * a,double * outi,double * outq,int N, int modu_index);
void AWGN(double *send, double *r, double sigma, int totallength);
void demodule(double *symbol_i, double *symbol_q, int symbol_len,double* out,double Kf,int modu_index);
//void DePuncture(double *received_punced_source, double *flow_for_decode, int source_length);
//void DeRepetition(double *received_punced_source, double *flow_for_decode, int source_length);
void TurboDecoding(double *flow_for_decode, int *flow_decoded,int flow_length);
void TurboDecodingfixpoint(double *flow_for_decode, int *flow_decoded,int flow_length);
void TurboCodingRelease();
void rate_match(int* input,int in_len,int*output,int out_len);
void de_rate_match(double*input,double*output,int in_len,int out_len);
