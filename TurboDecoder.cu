// TurboDecoder : Defines the entry point for the console application.
#include "helper_cuda.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <iostream>
using namespace std;


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.7E-308 //MINDOUBLE
#define RNMX (1.0-EPS)
#define INIFINITY  1E+10

#define MIN 1E-300
#define L_TOTAL 6144 // if u want to use block interleave,L_TOTAL must = x^2
#define MAXITER 10
#define	FRAME_NUM 10
#define AlphaBetaTHREAD_NUM 4

#define THREAD_NUM 768
#define BLOCK_NUM 16

typedef unsigned char BYTE;
typedef int INT;
typedef unsigned int UINT;
typedef int BOOL;

UINT m_Inter_table[L_TOTAL];
UINT LTE_Inter_table[L_TOTAL][3] =
	{40,3,10,	48,7,12,    56,19,42,	64,7,16,    72,7,18,    80,11,20,   
    88,5,22,    96,11,24,   104,7,26,   112,41,84,  120,103,90, 128,15,32,  
    136,9,34,   144,17,108, 152,9,38,	160,21,120, 168,101,84, 176,21,44,  
    184,57,46,  192,23,48,  200,13,50,  208,27,52,  216,11,36,  224,27,56,  
    232,85,58,  240,29,60,  248,33,62,  256,15,32,  264,17,198, 272,33,68,  
    280,103,210,288,19,36,  296,19,74,  304,37,76,  312,19,78,  320,21,120, 
    328,21,82,  336,115,84, 344,193,86, 352,21,44,  360,133,90, 368,81,46,  
    376,45,94,  384,23,48,  392,243,98, 400,151,40, 408,155,102,416,25,52,  
    424,51,106, 432,47,72,  440,91,110, 448,29,168, 456,29,114, 464,247,58, 
    472,29,118, 480,89,180, 488,91,122, 496,157,62, 504,55,84,  512,31,64,  
    528,17,66,  544,35,68,  560,227,420,576,65,96,  592,19,74,  608,37,76,  
    624,41,234, 640,39,80,  656,185,82, 672,43,252, 688,21,86,  704,155,44, 
    720,79,120, 736,139,92, 752,23,94,  768,217,48, 784,25,98,  800,17,80,  
    816,127,102,832,25,52,  848,239,106,864,17,48,  880,137,110,896,215,112,
    912,29,114, 928,15,58,  944,147,118,960,29,60,  976,59,122, 992,65,124, 
    1008,55,84, 1024,31,64, 1056,17,66, 1088,171,204,1120,67,140,1152,35,72,
    1184,19,74, 1216,39,76, 1248,19,78, 1280,199,240,1312,21,82,1344,211,252,
    1376,21,86, 1408,43,88, 1440,149,60,1472,45,92, 1504,49,846,1536,71,48, 
    1568,13,28, 1600,17,80, 1632,25,102,1664,183,104,1696,55,954,1728,127,96,
    1760,27,110,1792,29,112,1824,29,114,1856,57,116,1888,45,354,1920,31,120,
    1952,59,610,1984,185,124,2016,113,420,2048,31,64,2112,17,66,2176,171,136,
    2240,209,420,2304,253,216,2368,367,444,2432,265,456,2496,181,468,2560,39,80,
    2624,27,164,2688,127,504,2752,143,172,2816,43,88,2880,29,300,2944,45,92,
    3008,157,188,3072,47,96,3136,13,28,3200,111,240,3264,443,204,3328,51,104,
    3392,51,212,3456,451,192,3520,257,220,3584,57,336,3648,313,228,3712,271,232,
    3776,179,236,3840,331,120,3904,363,244,3968,375,248,4032,127,168,4096,31,64,
    4160,33,130,4224,43,264,4288,33,134,4352,477,408,4416,35,138,4480,233,280,
    4544,357,142,4608,337,480,4672,37,146,4736,71,444,4800,71,120,4864,37,152,
    4928,39,462,4992,127,234,5056,39,158,5120,39,80,5184,31,96,5248,113,902,
    5312,41,166,5376,251,336,5440,43,170,5504,21,86,5568,43,174,5632,45,176,
    5696,45,178,5760,161,120,5824,89,182,5888,323,184,5952,47,186,6016,23,94,
    6080,47,190,6144,263,480};

//typedef enum __bool { false = 0, true = 1, } bool;

long idum2;
long idum;
long iy;
long iv[NTAB];	
unsigned memory;

/*
Long period (? 2 \Theta 10 18 ) random number generator of L'Ecuyer with Bays­Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.
--
*/


double ran2()
{
	int j;
	long k;
	double temp;
	
	
	k=(idum)/IQ1;
	idum=IA1*(idum-k*IQ1)-k*IR1;  // Compute idum=(IA1*idum) % IM1 without overflows by Schrage's method.
	if (idum < 0)
		idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;	// Compute idum2=(IA2*idum) % IM2 likewise.
	if (idum2 < 0)
		idum2 += IM2;
	j = iy/NDIV;
	iy=iv[j]-idum2;
	//	iy=iv[j]-idum2; 		// Here idum is shuffled, idum and idum2 are combined to generate output.
	iv[j] = idum;
	if (iy < 1)
		iy += IMM1;
	if ((temp=AM*iy) > RNMX)
		return RNMX; 		// Because users don't expect endpoint values.
	else
		return temp;
}


void initRandom(long seed)
{
	idum2=123456789;
	idum=0;
	iy=0;
	
	if (seed != 0)
		idum = seed;
	else
		idum = 1;
	
	
    int j;
	for (j=NTAB+7;j>=0;j--) // Load the shuffle table (after 8 warm­ups).
	{
		long k=(idum)/IQ1;
		
		idum=IA1*(idum-k*IQ1)-k*IR1;
		if (idum < 0)
			idum += IM1;
		if (j < NTAB)
			iv[j] = idum;
	}
	iy=iv[0];
}


double doublerandom()
{
	double t = ran2();
	return t;
}

long longrandom(long range)
{
	double t;
	
	t = doublerandom();
	return((long)(t*(double)range));
}


bool boolrandom()
{
	double t=doublerandom();
	
	if (t>0.5)
		return true;
	else
		return false;
}
//time_t t;	
//time(&t);	
//init((long)t);
long seed = 1234421;





#define M	3	// register length,=tail length
#define NSTATE	8	// = M^2
#define L_ALL 3*L_TOTAL	// coded frame length
//#define DELTA 30	// SOVA window size. Make decision after 'delta' delay. Decide bit k when received bits
					// for bit (k+delta) are processed. Trace back from (k+delta) to k. 
// Rate 1/3 Turbo code
// The generator polynomials used are:
//	g0=[1 1 1]
//  g1=[1 0 1]
//	RSC encoder structure:
//
//			 +------------------------------------------>c1
//           |          g0(x)    
//           x_.--o-------------(+)<----------+
//           |    |             /|\           |
//			 |   \|/      ---    |     ---    |
// data --_.-o-->(+)--o->| D0|---o--->| D1|---o
//                    |   ---          ---    |
//                    |                       |
//                    +-------->(+)<----------+
//                      g1(x)    |
//								 +---------------------->c2
//

// NextOut[bk][current state]
static const char EnNextOut[2][NSTATE] = // check bit based on current and input bit
{	0,0,1,1,1,1,0,0,
	0,0,1,1,1,1,0,0
};
static const char NextOut[2][NSTATE] = // check bit based on current and input bit
{	-1,-1,1,1,1,1,-1,-1,
	-1,-1,1,1,1,1,-1,-1
};
// NextState[bk][current state]
static const BYTE NextState[2][NSTATE] = // next state based on current and input bit
{	0,4,5,1,2,6,7,3,
	4,0,1,5,6,2,3,7
};
// LastOut[bk][current state]
static const char LastOut[2][NSTATE] =	// trellis last check bit
{	-1,1,1,-1,-1,1,1,-1,
	-1,1,1,-1,-1,1,1,-1
};
// LastState[bk][current state]
static const BYTE LastState[2][NSTATE] =	// last state lead to current state by input bk
{	0,3,4,7,1,2,5,6,
	1,2,5,6,0,3,4,7
};
// TailBit[current state]
static const char TailBit[NSTATE] = // tail info bits when trellis is terminating
{	0,1,1,0,0,1,1,0
};








double gaussian(double variance)
{
	// static becuase we don't want to have it initialized each time we go in
	double returnvalue=0;
	double k;
	
	k = sqrt(variance/2.0);
	
	// add 24 uniform RV to obtain a simulation of normality
    int x;
	for (x=0;x<24;x++)
		returnvalue += doublerandom();
	
	return k*(returnvalue-0.5*24);

}




//////////////////////////////////////////////////////////////////////
// block interleave
// L_TOTAL must = x^2,otherwise,who knows?
//////////////////////////////////////////////////////////////////////
void init_Block_interleave_table()
{
	INT i,j;
	INT temp;

	temp = (INT)sqrt(L_TOTAL);
	for (i=0;i<temp;i++)
		for (j=0;j<temp;j++)
			m_Inter_table[i*temp+j] = j*temp+i;

	
}
void init_my_interleave_table() {
	INT i;
	for (i=0;i<L_TOTAL;i++)
		m_Inter_table[i] = L_TOTAL-1-i;
}

//////////////////////////////////////////////////////////////////////
// RSC endcoder
// mesg -- {0,1}
// parity -- {0,1}
// force==1,terminated --- for outer encoder
//////////////////////////////////////////////////////////////////////
void RSC_Encode(BYTE *mesg, BYTE *parity, unsigned int size, bool force)
{
	BYTE state,uk;
	unsigned x;
	
	state=0;
	for (x=0;x<size;x++)
	{
		// force the encoder to zero state at the end
		if (x>=size-M && force)
		{
			mesg[x] = TailBit[state];
		}
		
		// can't assume the bool type has an intrinsic value of 0 or 1
		// may differ from platform to platform
		uk = mesg[x] ? 1 : 0;
		
		// calculate output due to new mesg bit
		parity[x] = EnNextOut[uk][state];
		// calculate the new state
		state = NextState[uk][state];
	}
}


//////////////////////////////////////////////////////////////////////
// Turbo encoder
// msg -- {0,1}
// stream -- {0,1}
// puncture -- true to get 1/2 rate,NOT tested yet
//////////////////////////////////////////////////////////////////////
void encode(BYTE *msg, BYTE *stream, bool puncture)
{
	INT i;
	BYTE imsg[L_TOTAL];
	BYTE chkBuffer[2][L_TOTAL];
	// first encoder
	RSC_Encode(msg,chkBuffer[0],L_TOTAL,true);
	// interleave
	for (i=0;i<L_TOTAL;i++)
		imsg[i]=msg[m_Inter_table[i]];
	// second encoder
	RSC_Encode(imsg,chkBuffer[1],L_TOTAL,false);
	// punture
	for (i=0;i<L_TOTAL;i++)
	{
		if(puncture){
			stream[i*2]=msg[i];
			stream[i*2+1]=chkBuffer[i%2][i];
		}else{
			stream[i*3]=msg[i];
			stream[i*3+1]=chkBuffer[0][i];
			stream[i*3+2]=chkBuffer[1][i];
		}	
	}
}


__global__ void interLeave(double * src, double * des , unsigned int * interLeaveTable[3] ){
//__global__ void interLeave(double * src, double * des , unsigned int * interLeaveTable ){
    const int tid = threadIdx.x;
   //des[tid] = src[interLeaveTable[tid]];
    des[tid] = src[(interLeaveTable[tid][1]*tid+interLeaveTable[tid][2]*tid*tid)%interLeaveTable[tid][0]];
}

__global__ void deInterLeave(double * src, double * des , unsigned int * interLeaveTable[3] ){
//__global__ void deInterLeave(double * src, double * des , unsigned int * interLeaveTable ){
    const int tid = threadIdx.x;
    //des[interLeaveTable[tid]] = src[tid];
    des[(interLeaveTable[tid][1]*tid+interLeaveTable[tid][2]*tid*tid)%interLeaveTable[tid][0]] = src[tid];
}

__global__ void gammaAlpha(double * msg ,double * parity, double * L_a, double (*gamma)[8][8], BYTE (*lastState)[8],char (*lastOut)[8] ){
    const int tid = threadIdx.x;

    unsigned int s0,s2;
    for (s0=0;s0<NSTATE;s0++) {
		for (s2=0;s2<NSTATE;s2++)
			gamma[tid][s0][s2]=-INIFINITY;
		gamma[tid][s0][lastState[0][s0]]=-msg[tid]+parity[tid]*lastOut[0][s0]-log(1+exp(L_a[tid]));
		gamma[tid][s0][lastState[1][s0]]=msg[tid]+parity[tid]*lastOut[1][s0]+L_a[tid]-log(1+exp(L_a[tid]));
		//gamma[tid][s0][lastState[0][s0]]=0.5;
		//gamma[tid][s0][lastState[1][s0]]=-0.5;
    }
}

__global__ void gammaBeta(double * msg ,double * parity, double * L_a, double (*gamma)[8][8], BYTE (*nextState)[8], char (*nextOut)[8]){
    const int tid = threadIdx.x;

    unsigned int s0,s2;
    for (s0=0;s0<NSTATE;s0++) {
		for (s2=0;s2<NSTATE;s2++)
			gamma[tid][s0][s2]=-INIFINITY;
		gamma[tid][s0][nextState[0][s0]]=-msg[tid]+parity[tid]*nextOut[0][s0]-log(1+exp(L_a[tid]));
		gamma[tid][s0][nextState[1][s0]]=msg[tid]+parity[tid]*nextOut[1][s0]+L_a[tid]-log(1+exp(L_a[tid]));
		//gamma[tid][s0][nextState[0][s0]]=0.5;
		//gamma[tid][s0][nextState[1][s0]]=-0.5;
    }
}

__global__ void Alpha(double (*Alpha)[8], double (*gamma)[8][8], double *maxBranch) {
	const int tid = threadIdx.x;

	UINT k, s1, s2;
	double sum;

	if (tid == 0) {
		Alpha[0][0] = 0.0;
		for (s1=1;s1<NSTATE;s1++)
			Alpha[0][s1]=-INIFINITY;
	}
	else {
		for (s1=0;s1<NSTATE;s1++)
			Alpha[tid*(L_TOTAL/AlphaBetaTHREAD_NUM)][s1]=0;
	}

	//for (k=1; k<=L_TOTAL; k++) {
	for (k=tid*L_TOTAL/AlphaBetaTHREAD_NUM+1; k<(tid*L_TOTAL/AlphaBetaTHREAD_NUM+L_TOTAL/AlphaBetaTHREAD_NUM); k++) {
        for (s2=0;s2<NSTATE;s2++){
            sum = 0.0;
            for (s1=0;s1<NSTATE;s1++) {
                sum+=exp(gamma[k-1][s2][s1]+Alpha[k-1][s1]);
			}
            if (sum<MIN)
            //if (sum<=0.000000000000000000000000000001)
                Alpha[k][s2]=-INIFINITY;
            else
                Alpha[k][s2]=log(sum);
        }

		// normalization,prevent overflow
		maxBranch[k]=Alpha[k][0];
		for (s2=1;s2<NSTATE;s2++)
			if (Alpha[k][s2]>maxBranch[k])
				maxBranch[k]=Alpha[k][s2];

		for (s2=0;s2<NSTATE;s2++)
			Alpha[k][s2]=Alpha[k][s2]-maxBranch[k];
	}

}

__global__ void Beta(double (*Beta)[8], double (*gamma)[8][8], bool index, double* maxBranch) {
	const int tid = threadIdx.x;

	UINT k, s1, s2;
	double sum;

	if (tid == (AlphaBetaTHREAD_NUM-1)) {
		if (index){// true -- terminated,false -- open
        Beta[L_TOTAL][0]=0.0;
        for (s2=1;s2<NSTATE;s2++)
            Beta[L_TOTAL][s2]=-INIFINITY;
		}
		else 
			for (s2=0;s2<NSTATE;s2++)
				Beta[L_TOTAL][s2]=0;
	}
	else {
		for (s2=0; s2<NSTATE; s2++)
			Beta[(tid+1)*L_TOTAL/AlphaBetaTHREAD_NUM][s2]=0;
	}

    for (k=(tid+1)*L_TOTAL/AlphaBetaTHREAD_NUM-1;k>(tid*L_TOTAL/AlphaBetaTHREAD_NUM);k--) {
   // for (k=L_TOTAL-1;k>0;k--) {

        for (s1=0;s1<NSTATE;s1++) {
            sum = 0.0;
            for (s2=0;s2<NSTATE;s2++) 
                sum += exp(gamma[k][s1][s2] + Beta[k+1][s2]);
            if (sum<MIN)
            //if (sum<=0.000000000000000000000000000001)
                Beta[k][s1] = -INIFINITY;
            else 
                Beta[k][s1] = log(sum);
        }

		// normalization,prevent overflow
		for (s2=0;s2<NSTATE;s2++)
			Beta[k][s2]=Beta[k][s2]-maxBranch[k];
	}
}

void computeAlpha(double (*AlphaHost)[8], double (*gamma)[8][8], double *maxBranch) {
    // initialize Alpha & Beta
    AlphaHost[0][0]=0;
	UINT s1,k,s2;
	double sum;
    for (s1=1;s1<NSTATE;s1++)
        AlphaHost[0][s1]=-INIFINITY;

    for (k=1;k<=L_TOTAL;k++){

        for (s2=0;s2<NSTATE;s2++){
            sum = 0;
            for (s1=0;s1<NSTATE;s1++) {
                sum+=exp(gamma[k-1][s2][s1]+AlphaHost[k-1][s1]);
			}
            if (sum<MIN)
            //if (sum<=0.000000000000000000000000000001)
                AlphaHost[k][s2]=-INIFINITY;
            else
                AlphaHost[k][s2]=log(sum);
        }

		// normalization,prevent overflow
		maxBranch[k]=AlphaHost[k][0];
		for (s2=1;s2<NSTATE;s2++)
			if (AlphaHost[k][s2]>maxBranch[k])
				maxBranch[k]=AlphaHost[k][s2];

		for (s2=0;s2<NSTATE;s2++)
			AlphaHost[k][s2]=AlphaHost[k][s2]-maxBranch[k];
    }
}

void computeBeta(double (*BetaHost)[8], double (*gamma)[8][8], bool index, double * maxBranch){
    // initialize Beta
	UINT s1,k,s2;
	double sum;
    if (index){// true -- terminated,false -- open
        BetaHost[L_TOTAL][0]=0;
        for (s2=1;s2<NSTATE;s2++)
            BetaHost[L_TOTAL][s2]=-INIFINITY;
    }
    else 
        for (s2=0;s2<NSTATE;s2++)
            BetaHost[L_TOTAL][s2]=0;

    for (k=L_TOTAL-1;k>0;k--) {

        for (s1=0;s1<NSTATE;s1++) {
            sum = 0.0;
            for (s2=0;s2<NSTATE;s2++) 
                sum += exp(gamma[k][s1][s2] + BetaHost[k+1][s2]);
            if (sum<MIN)
            //if (sum<=0.000000000000000000000000000001)
                BetaHost[k][s1] = -INIFINITY;
            else 
                BetaHost[k][s1] = log(sum);
        }

		// normalization,prevent overflow
		for (s2=0;s2<NSTATE;s2++)
			BetaHost[k][s2]=BetaHost[k][s2]-maxBranch[k];
    }
}

__global__ void normalizationAlphaAndBeta(double (*Alpha)[8], double (*Beta)[8]) {
    unsigned int tid = threadIdx.x+1; 
    double max_branch;
    max_branch = Alpha[tid][0];
	UINT s2;
    for (s2=1;s2<NSTATE;s2++)
        if (Alpha[tid][s2]>max_branch)
            max_branch = Alpha[tid][s2];

    for (s2=0;s2<NSTATE;s2++) {
        Alpha[tid][s2] = Alpha[tid][s2] - max_branch;

        if (tid != L_TOTAL) 
            Beta[tid][s2] = Beta[tid][s2] - max_branch;
    }

}

__global__ void LLRS(double * msg, double * parity, double * L_a, double (*Alpha)[8], double (*Beta)[8], double * L_all,BYTE (*lastState)[8], char (*lastOut)[8]) {
    unsigned int tid = threadIdx.x; 
    UINT s2;
	double sum0 = 0.0, sum1 = 0.0;
    for (s2=0;s2<NSTATE;s2++) {
        //gamma[LastState[0][s2]]=-msg[tid]+parity[tid]*LastOut[0][s2]-log(1+exp(L_a[tid]));
        //gamma[LastState[1][s2]]=msg[tid]+parity[tid]*LastOut[1][s2]+L_a[tid]-log(1+exp(L_a[tid]));
        //sum0+=exp(gamma[LastState[0][s2]]+Alpha[tid][LastState[0][s2]]+Beta[tid+1][s2]);
        //sum1+=exp(gamma[LastState[1][s2]]+Alpha[tid][LastState[1][s2]]+Beta[tid+1][s2]);
        double gamma0=-msg[tid]+parity[tid]*lastOut[0][s2]-log(1+exp(L_a[tid]));
        double gamma1=msg[tid]+parity[tid]*lastOut[1][s2]+L_a[tid]-log(1+exp(L_a[tid]));
        sum0+=exp(gamma0+Alpha[tid][lastState[0][s2]]+Beta[tid+1][s2]);
        sum1+=exp(gamma1+Alpha[tid][lastState[1][s2]]+Beta[tid+1][s2]);
    }
    //L_all[tid]=log(sum1)-log(sum0);
    L_all[tid]=log(sum1)-log(sum0);
}

__global__ void extrinsicInformation(double * L_all, double * msg, double * L_a, double * L_e) {
    unsigned int tid = threadIdx.x;
    L_e[tid] = L_all[tid] - 2*msg[tid] - L_a[tid];
}

__global__ void demultiplex(double * stream, double * msg, double * parity0, double * parity1) {
    unsigned int tid = threadIdx.x;
    //if (puncture){// punctured rate=1/2
    //    msg[tid]=stream[2*tid];
    //    parity[tid%2][tid]=stream[tid*2+1];
    //}
    //else {// unpunctured rate=1/3
    //    msg[tid]=stream[3*tid];
    //    parity0[tid]=stream[3*tid+1];
    //    parity1[tid]=stream[3*tid+2];
    //}
        msg[tid]=stream[3*tid];
        parity0[tid]=stream[3*tid+1];
        parity1[tid]=stream[3*tid+2];
}

__global__ void initializeExtrinsicInformation(double * L_e) {
    unsigned int tid = threadIdx.x;
    L_e[tid] = 0;
    
}

__global__ void exestimateInformationBits(double * L_all, BYTE * msghat, UINT * interLeaveTable[3]) {
//__global__ void exestimateInformationBits(double * L_all, BYTE * msghat, UINT * m_Inter_table) {
    unsigned int tid = threadIdx.x;
    if(L_all[tid]>0)
		msghat[(interLeaveTable[tid][1]*tid+interLeaveTable[tid][2]*tid*tid)%interLeaveTable[tid][0]]=1;
        //msghat[m_Inter_table[tid]]=1;
    else
		msghat[(interLeaveTable[tid][1]*tid+interLeaveTable[tid][2]*tid*tid)%interLeaveTable[tid][0]]=0;
        //msghat[m_Inter_table[tid]]=0;
}

void countErrors(BYTE *m, BYTE * mhat, UINT * bitsError, UINT * frameError, UINT iter) {

	bool f_err = false;
	for (int i=0; i<(L_TOTAL-M);i++) {
		if (m[i] != mhat[i]) {
			bitsError[iter] = bitsError[iter]+1;
			f_err = true;
		}
	}

	if (f_err) 
		frameError[iter] = frameError[iter]+1;
}


int main(int argc, char* argv[])
{
    initRandom(seed);

	BYTE * m;
	BYTE * x;
	double * y;
	BYTE * mhat;

	int frame;
	UINT bits_all,bits_err[MAXITER],frame_err[MAXITER];
	double Ber,Fer;
	double Eb_No_dB,No;
	//bool f_err;
	//FILE * fp;
	int i;

	m = new BYTE[L_TOTAL];
	x = new BYTE[L_ALL];
	y = new double[L_ALL];
	mhat = new BYTE[L_TOTAL];

	//init_Block_interleave_table();	// block interleave
	//init_my_interleave_table();


    
    findCudaDevice(argc, (const char **)argv);

	BYTE (*LastStateDevice)[8];
	BYTE (*NextStateDevice)[8];
	char (*LastOutDevice)[8];
	char (*NextOutDevice)[8];
	double * yDevice;
	double * msgDevice;
	double * imsgDevice;
	BYTE * mhatDevice;
	double * parity0Device;
	double * parity1Device;
	//UINT * tableDevice;
	UINT * LTETableDevice[3];
	double * L_eDevice;
	double * L_aDevice;
	double * L_allDevice;
	double (*gammaAlphaDevice)[8][8];
	double (*gammaBetaDevice)[8][8];
	double (*AlphaDevice)[8];
	double (*BetaDevice)[8];
	double *maxBranchDevice;

    cudaMalloc((void **)&LastStateDevice, 2*8*sizeof(BYTE));
    cudaMalloc((void **)&NextStateDevice, 2*8*sizeof(BYTE));
    cudaMalloc((void **)&LastOutDevice, 2*8*sizeof(char));
    cudaMalloc((void **)&NextOutDevice, 2*8*sizeof(char));

    cudaMalloc((void **)&yDevice, L_ALL*sizeof(double));
    cudaMalloc((void **)&msgDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&imsgDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&mhatDevice, L_TOTAL*sizeof(BYTE));
    cudaMalloc((void **)&parity0Device, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&parity1Device, L_TOTAL*sizeof(double));
    //cudaMalloc((void **)&tableDevice, L_TOTAL*sizeof(unsigned int));
    cudaMalloc((void **)&LTETableDevice, L_TOTAL*3*sizeof(unsigned int));
    cudaMalloc((void **)&L_eDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&L_aDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&L_allDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&gammaAlphaDevice, L_TOTAL*sizeof(double)*8*8);
    cudaMalloc((void **)&gammaBetaDevice, L_TOTAL*sizeof(double)*8*8);
    cudaMalloc((void **)&AlphaDevice, (L_TOTAL+1)*sizeof(double)*8);
    cudaMalloc((void **)&BetaDevice, (L_TOTAL+1)*sizeof(double)*8);
    cudaMalloc((void **)&maxBranchDevice, (L_TOTAL+1)*sizeof(double));

	//For Debug
	//double L_aHost[L_TOTAL];
	//double L_aHost1[L_TOTAL];
	//double L_allHost[L_TOTAL];

    double gammaAlphaHost[L_TOTAL][8][8];
    double gammaBetaHost[L_TOTAL][8][8];
    double AlphaHost[L_TOTAL+1][8];
    double BetaHost[L_TOTAL+1][8];

	double max_branch[L_TOTAL+1];

    cudaMemcpy(LastStateDevice,LastState,sizeof(BYTE)*2*8, cudaMemcpyHostToDevice);
    cudaMemcpy(NextStateDevice,NextState,sizeof(BYTE)*2*8, cudaMemcpyHostToDevice);
    cudaMemcpy(LastOutDevice,LastOut,sizeof(char)*2*8, cudaMemcpyHostToDevice);
    cudaMemcpy(NextOutDevice,NextOut,sizeof(char)*2*8, cudaMemcpyHostToDevice);

    //cudaMemcpy(tableDevice,m_Inter_table,sizeof(unsigned int)*L_TOTAL, cudaMemcpyHostToDevice);
    cudaMemcpy(LTETableDevice,LTE_Inter_table,sizeof(unsigned int)*L_TOTAL*3, cudaMemcpyHostToDevice);

	for (Eb_No_dB= 0.0;Eb_No_dB<1.0;Eb_No_dB+=1){

	//Eb_No_dB = 0.0;
		No = 1/pow(10.0,Eb_No_dB/10.0);
		bits_all = 0;
		for (i =0; i<MAXITER;i++) {
			bits_err[i]=0;
			frame_err[i]=0;
		}

		for (frame = 0; frame<FRAME_NUM; frame++, bits_all += (L_TOTAL-M)) {

			// Generate random information bits
			for (i=0;i<L_TOTAL;i++)
				if (boolrandom())
					m[i]=1;
				else
					m[i]=0;
			// encoder
			encode(m,x,false);
			// add noise
			for (i=0;i<L_ALL;i++)
				if (x[i])
					y[i]=1.0+gaussian(No/2);
				else
					y[i]=-1.0+gaussian(No/2);

			cudaMemcpy(yDevice,y,sizeof(double)*L_ALL, cudaMemcpyHostToDevice);

			demultiplex<<<1,L_TOTAL>>>(yDevice, msgDevice, parity0Device, parity1Device); 
			//interLeave<<<1,L_TOTAL>>>(msgDevice, imsgDevice, tableDevice);
			interLeave<<<1,L_TOTAL>>>(msgDevice, imsgDevice, LTETableDevice);
			initializeExtrinsicInformation<<<1,L_TOTAL>>>(L_eDevice);

			for (int iter = 0; iter<MAXITER; iter++) {
				
			//	deInterLeave<<<1,L_TOTAL>>>(L_eDevice, L_aDevice, tableDevice);
				deInterLeave<<<1,L_TOTAL>>>(L_eDevice, L_aDevice, LTETableDevice);

				gammaAlpha<<<1,L_TOTAL>>>(msgDevice , parity0Device,  L_aDevice,  gammaAlphaDevice,LastStateDevice, LastOutDevice);
				gammaBeta<<<1,L_TOTAL>>>(msgDevice , parity0Device,  L_aDevice,  gammaBetaDevice, NextStateDevice, NextOutDevice);
				Alpha<<<1,AlphaBetaTHREAD_NUM>>>(AlphaDevice, gammaAlphaDevice, maxBranchDevice);
				Beta<<<1,AlphaBetaTHREAD_NUM>>>(BetaDevice, gammaBetaDevice,true, maxBranchDevice);
				cudaMemcpy(AlphaHost, AlphaDevice, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyDeviceToHost);
				cudaMemcpy(gammaAlphaHost, gammaAlphaDevice, sizeof(double)*L_TOTAL*8*8, cudaMemcpyDeviceToHost);
				//cudaMemcpy(gammaBetaHost, gammaBetaDevice, sizeof(double)*L_TOTAL*8*8, cudaMemcpyDeviceToHost);

				//computeAlpha(AlphaHost, gammaAlphaHost, max_branch);
				//computeBeta(BetaHost, gammaBetaHost, true,max_branch);
				//cudaMemcpy(AlphaDevice, AlphaHost, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyHostToDevice);
				//cudaMemcpy(BetaDevice, BetaHost, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyHostToDevice);
				normalizationAlphaAndBeta<<<1,L_TOTAL>>>(AlphaDevice, BetaDevice);
				cudaMemcpy(AlphaHost, AlphaDevice, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyDeviceToHost);

				LLRS<<<1,L_TOTAL>>>(msgDevice, parity0Device, L_aDevice, AlphaDevice, BetaDevice, L_allDevice,LastStateDevice, LastOutDevice);

				extrinsicInformation<<<1, L_TOTAL>>>(L_allDevice, msgDevice, L_aDevice, L_eDevice);
				//if (iter >= 3) {
				///debug
				//exestimateInformationBits<<<1,L_TOTAL>>>(L_allDevice, mhatDevice, tableDevice); 

				//cudaMemcpy(mhat, mhatDevice, sizeof(BYTE)*L_TOTAL, cudaMemcpyDeviceToHost);
				//countErrors(m, mhat, bits_err, frame_err, iter);
				//debug
				//cudaMemcpy(L_aHost1, L_aDevice, sizeof(double)*L_TOTAL, cudaMemcpyDeviceToHost);
				//}

				//interLeave<<<1, L_TOTAL>>>(L_eDevice, L_aDevice, tableDevice);
				interLeave<<<1, L_TOTAL>>>(L_eDevice, L_aDevice, LTETableDevice);

				gammaAlpha<<<1,L_TOTAL>>>(imsgDevice , parity1Device,  L_aDevice,  gammaAlphaDevice, LastStateDevice, LastOutDevice);
				gammaBeta<<<1,L_TOTAL>>>(imsgDevice , parity1Device,  L_aDevice,  gammaBetaDevice, NextStateDevice, NextOutDevice);
				Alpha<<<1,AlphaBetaTHREAD_NUM>>>(AlphaDevice, gammaAlphaDevice, maxBranchDevice);
				Beta<<<1,AlphaBetaTHREAD_NUM>>>(BetaDevice, gammaBetaDevice,false, maxBranchDevice);
				cudaMemcpy(AlphaHost, AlphaDevice, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyDeviceToHost);
				//cudaMemcpy(gammaAlphaHost, gammaAlphaDevice, sizeof(double)*L_TOTAL*8*8, cudaMemcpyDeviceToHost);
				//cudaMemcpy(gammaBetaHost, gammaBetaDevice, sizeof(double)*L_TOTAL*8*8, cudaMemcpyDeviceToHost);

				//computeAlpha(AlphaHost, gammaAlphaHost, max_branch);
				//Beta(BetaHost, gammaBetaHost, false, max_branch);
				//cudaMemcpy(AlphaDevice, AlphaHost, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyHostToDevice);
				//cudaMemcpy(BetaDevice, BetaHost, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyHostToDevice);
				normalizationAlphaAndBeta<<<1,L_TOTAL>>>(AlphaDevice, BetaDevice);
				cudaMemcpy(AlphaHost, AlphaDevice, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyDeviceToHost);

				LLRS<<<1,L_TOTAL>>>(imsgDevice, parity1Device, L_aDevice, AlphaDevice, BetaDevice, L_allDevice, LastStateDevice,LastOutDevice);

				extrinsicInformation<<<1, L_TOTAL>>>(L_allDevice, imsgDevice, L_aDevice, L_eDevice);

				exestimateInformationBits<<<1,L_TOTAL>>>(L_allDevice, mhatDevice, LTETableDevice); 
				//exestimateInformationBits<<<1,L_TOTAL>>>(L_allDevice, mhatDevice, tableDevice); 

				cudaMemcpy(mhat, mhatDevice, sizeof(BYTE)*L_TOTAL, cudaMemcpyDeviceToHost);
				countErrors(m, mhat, bits_err, frame_err, iter);

				//debug
				//cudaMemcpy(L_aHost, L_aDevice, sizeof(double)*L_TOTAL, cudaMemcpyDeviceToHost);
			}
			// estimate information bits
			//exestimateInformationBits<<<1,L_TOTAL>>>(L_allDevice, mhatDevice, tableDevice); 

			//cudaMemcpy(mhat, mhatDevice, sizeof(BYTE)*L_TOTAL, cudaMemcpyDeviceToHost);
			// count errors
			//UINT bits_err = 0;
			//for (i=0;i<L_TOTAL-M;i++) {
			//	if (mhat[i]!=m[i]) {
			//		bits_err++;
			//	}
			//}
			//cout<<"bits_err: "<<bits_err<<endl;
		}

		printf("-------------------------\n");
		printf("Eb/No=%fdB:\n",Eb_No_dB);
		printf("-------------------------\n");
		//fprintf(fp,"-------------------------\n");
		//fprintf(fp,"Eb/No=%fdB:\n",Eb_No_dB);
		//fprintf(fp,"-------------------------\n");

		for (i=0;i<MAXITER;i++) {
			Ber=(double)bits_err[i]/(double)bits_all;
			Fer=(double)frame_err[i]/(double)FRAME_NUM;
			printf("Iteration:%d\n",i+1);
			printf("---Ber=%f\n---Fer=%f\n",Ber,Fer);
			//fprintf(fp,"Iteration:%d\n",i);
			//fprintf(fp,"---Ber=%f\n---Fer=%f\n",Ber,Fer);
		}
	}

	delete m;
	delete x;
	delete y;
	delete mhat;

	cudaFree(LastStateDevice);
	cudaFree(NextStateDevice);
	cudaFree(LastOutDevice);
	cudaFree(NextOutDevice);
	cudaFree(yDevice);
	cudaFree(msgDevice);
	cudaFree(imsgDevice);
	cudaFree(mhatDevice);
	cudaFree(parity0Device);
	cudaFree(parity1Device);
//	cudaFree(tableDevice);
	cudaFree(LTETableDevice);
	cudaFree(L_eDevice);
	cudaFree(L_aDevice);
	cudaFree(L_allDevice);
	cudaFree(gammaAlphaDevice);
	cudaFree(gammaBetaDevice);
	cudaFree(AlphaDevice);
	cudaFree(BetaDevice);

	return 0;
}
