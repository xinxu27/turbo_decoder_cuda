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




#define L_TOTAL 128	// if u want to use block interleave,L_TOTAL must = x^2
#define M	2	// register length,=tail length
#define NSTATE	4	// = M^2
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
typedef unsigned char BYTE;
typedef int INT;
typedef unsigned int UINT;
typedef int BOOL;

// NextOut[bk][current state]
static const char EnNextOut[2][NSTATE] = // check bit based on current and input bit
{	0,0,1,1,
	1,1,0,0
};
static const char NextOut[2][NSTATE] = // check bit based on current and input bit
{	-1,-1,1,1,
	1,1,-1,-1
};
// NextState[bk][current state]
static const BYTE NextState[2][NSTATE] = // next state based on current and input bit
{	0,2,3,1,
	2,0,1,3
};
// LastOut[bk][current state]
static const char LastOut[2][NSTATE] =	// trellis last check bit
{	-1,1,-1,1,
	1,-1,1,-1
};
// LastState[bk][current state]
static const BYTE LastState[2][NSTATE] =	// last state lead to current state by input bk
{	0,3,1,2,
	1,2,0,3
};
// TailBit[current state]
static const char TailBit[NSTATE] = // tail info bits when trellis is terminating
{	0,1,1,0
};

//#define MAXITER 5
//#define	FRAME_NUM 10






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


__global__ void interLeave(double * src, double * des , unsingned int * interLeaveTable ){
    const int tid = threadIdx.x;
    des[tid] = src(interLeaveTable[tid]);
}

__global__ void deInterLeave(double * src, double * des , unsingned int * interLeaveTable ){
    const int tid = threadIdx.x;
    des(interLeaveTable[tid]) = src[i];
}

__global__ void gammaAlpha(double * msg ,double * parity, double * L_a, double (*gamma)[4][4]){
    const int tid = threadIdx.x;

    unsigned int s0, s1,s2;
    for (s0=0;s0<NSTATE;s0++) {
        for (s1=0;s1<NSTATE;s1++) {
            for (s2=0;s2<NSTATE;s2++)
                gammaAlpha[tid][s0][s2]=-INIFINITY;
            gammaAlpha[tid][s0][LastStateDevice[0][s1]]=-msg[tid]+parity[tid]*LastOutDevice[0][s1]-log(1+exp(L_a[tid]));
            gammaAlpha[tid][s0][LastStateDevice[1][s1]]=msg[tid]+parity[tid]*LastOutDevice[1][s1]+L_a[tid]-log(1+exp(L_a[tid]));
        }
    }
}

__global__ void gammaBeta(double * msg ,double * parity, double * L_a, double (*gamma)[4][4]){
    const int tid = threadIdx.x;

    unsigned int s0, s1,s2;
    for (s0=0;s0<NSTATE;s0++) {
        for (s1=0;s1<NSTATE;s1++) {
            for (s2=0;s2<NSTATE;s2++)
                gammaBeta[tid][s0][s2]=-INIFINITY;
            gammaBeta[tid][s0][NextStateDevice[0][s1]]=-msg[tid]+parity[tid]*NextOutDevice[0][s1]-log(1+exp(L_a[tid]));
            gammaBeta[tid][s0][NextStateDevice[1][s1]]=msg[tid]+parity[tid]*NextOutDevice[1][s1]+L_a[tid]-log(1+exp(L_a[tid]));
        }
    }
}

void Alpha(double (*AlphaHost)[4], double * gamma[4][4]) {
    // initialize Alpha & Beta
    AlphaHost[0][0]=0;
    for (s1=1;s1<nstate;s1++)
        AlphaHost[0][s1]=-INIFINITY;

    for (k=1,k<=L_TOTAL;k++){

        for (s2=0;s2<nstate;s2++){
            sum = 0;
            for (s1=0;s1<nstate;s1++)
                sum+=exp(gamma[k-1][s2][s1]+AlphaHost[k-1][s1]);
            if (sum<1E-300)
                AlphaHost[k][s2]=-INIFINITY;
            else
                AlphaHost[k][s2]=log(sum);
        }
    }
}

void Beta(double (*BetaHost)[4], double (*gamma)[4][4]){
    // initialize Beta
    if (index){// true -- terminated,false -- open
        BetaHost[L_TOTAL][0]=0;
        for (s2=1;s2<nstate;s2++)
            BetaHost[L_TOTAL][s2]=-INIFINITY;
    }
    else 
        for (s2=0;s2<nstate;s2++)
            BetaHost[L_TOTAL][s2]=0;

    for (k=L_TOTAL-1;k>0;k--) {

        for (s1=0;s1<nstate;s1++) {
            sum = 0.0;
            for (s2=0;s2<nstate;s2++) 
                sum += exp(gamma[k][s1][s2] + BetaHost[k+1][s2]);
            if (sum<1E-300)
                BetaHost[k][s1] = -INIFINITY;
            else 
                BetaHost[k][s1] = log(sum);
        }
    }
}

__global__ void normalizationAlphaAndBeta(double (*Alpha)[4], double (*Beta)[4]) {
    unsigned int tid = threadIdx.x+1; 
    double max_branch[L_TOTAL];
    max_branch[tid] = Alpha[k][0];
    for (s2=1;s2<NSTATE;s2++)
        if (Alpha[tid][s2]>max_branch[tid])
            max_branch[tid] = Alpha[tid][s2];

    for (s2=0;s2<NSTATE;s2++) {
        Alpha[tid][s2] = Alpha[tid][s2] - max_branch[tid];

        if (tid != L_TOTAL) 
            Beta[tid][s2] = Beta[tid][s2] - max_branch[tid];
    }

}

__global__ void LLRS(double * msg, double * parity, double * L_a, double (*Alpha)[4], double (*Beta)[4], double * L_all) {
    unsigned int tid = threadIdx.x; 
    UINT sum0=0;sum1=0;
    for (s2=0;s2<NSTATE;s2++) {
        //gamma[LastState[0][s2]]=-msg[tid]+parity[tid]*LastOut[0][s2]-log(1+exp(L_a[tid]));
        //gamma[LastState[1][s2]]=msg[tid]+parity[tid]*LastOut[1][s2]+L_a[tid]-log(1+exp(L_a[tid]));
        //sum0+=exp(gamma[LastState[0][s2]]+Alpha[tid][LastState[0][s2]]+Beta[tid+1][s2]);
        //sum1+=exp(gamma[LastState[1][s2]]+Alpha[tid][LastState[1][s2]]+Beta[tid+1][s2]);
        double gamma0=-msg[tid]+parity[tid]*LastOutDevice[0][s2]-log(1+exp(L_a[tid]));
        double gamma1=msg[tid]+parity[tid]*LastOutDevice[1][s2]+L_a[tid]-log(1+exp(L_a[tid]));
        sum0+=exp(gamma0+Alpha[tid][LastStateDevice[0][s2]]+Beta[tid+1][s2]);
        sum1+=exp(gamma1+Alpha[tid][LastStateDevice[1][s2]]+Beta[tid+1][s2]);
    }
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

__global__ initializeExtrinsicInformation(double * L_e) {
    unsigned int tid = threadIdx.x;
    L_e[tid] = 0;
    
}

__global__ exestimateInformationBits(double * L_all, BYTE * msghat, UINT * m_Inter_table) {
    unsigned int tid = threadIdx.x;
    if(L_all[tid]>0)
        msghat[m_Inter_table[tid]]=1;
    else
        msghat[m_Inter_table[tid]]=0;
}


int main(int argc, char* argv[])
{
    initRandom(seed);

	BYTE * m;
	BYTE * x;
	double * y;
	BYTE * mhat;

	//int frame;
	//int bits_all,bits_err[MAXITER],frame_err[MAXITER];
	//double Ber,Fer;
	double Eb_No_dB,No;
	//bool f_err;
	//FILE * fp;
	int i;

	m = new BYTE[L_TOTAL];
	x = new BYTE[L_ALL];
	y = new double[L_ALL];
	mhat = new BYTE[L_TOTAL];

    UINT m_Inter_table[L_TOTAL];
	init_Block_interleave_table();	// block interleave


	Eb_No_dB=-3.0;
    No = 1/pow(10.0,Eb_No_dB/10.0);

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
    
    findCudaDevice(argc, (const char *)argv);

    cudaMalloc((void **)&LastStateDevice, 2*4*sizeof(BYTE));
    cudaMalloc((void **)&NextStateDevice, 2*4*sizeof(BYTE));
    cudaMalloc((void **)&LastOutDevice, 2*4*sizeof(char));
    cudaMalloc((void **)&NextOutDevice, 2*4*sizeof(char));

    cudaMalloc((void **)&yDevice, L_ALL*sizeof(double));
    cudaMalloc((void **)&msgDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&mhatDevice, L_TOTAL*sizeof(BYTE));
    cudaMalloc((void **)&parity0Device, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&parity1Device, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&tableDevice, L_TOTAL*sizeof(unsigned int));
    cudaMalloc((void **)&imsgDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&L_eDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&L_aDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&L_allDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&L_allDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&L_eDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&gammaAlphaDevice[4][4], L_TOTAL*sizeof(double)*4*4);
    cudaMalloc((void **)&gammaBetaDevice[4][4], L_TOTAL*sizeof(double)*4*4);
    cudaMalloc((void **)&AlphaDevice[4], L_TOTAL*sizeof(double)*4);
    cudaMalloc((void **)&BetaDevice[4], L_TOTAL*sizeof(double)*4);
    double gammaAlphaHost[L_TOTAL][4][4];
    double gammaBetaHost[L_TOTAL][4][4];
    double AlphaHost[L_TOTAL][4];
    double BetaHost[L_TOTAL][4];


    cudaMemcpy(LastStateDevice,LastState,sizeof(BYTE)*2*4, cudaMemcpyHostToDevice);
    cudaMemcpy(NextStateDevice,NextState,sizeof(BYTE)*2*4, cudaMemcpyHostToDevice);
    cudaMemcpy(LastOutDevice,LastOut,sizeof(char)*2*4, cudaMemcpyHostToDevice);
    cudaMemcpy(NextOutDevice,NextOut,sizeof(char)*2*4, cudaMemcpyHostToDevice);

    cudaMemcpy(tableDevice,m_Inter_table,sizeof(unsigned int)*L_TOTAL, cudaMemcpyHostToDevice);
    cudaMemcpy(yDevice,y,sizeof(double)*L_ALL, cudaMemcpyHostToDevice);

    demultiplex<<<1,L_TOTAL>>>(yDevice, msgDevice, parity0Device, parity1Device); 
    interLeave<<<1,L_TOTAL>>>(msgDevice, imsgDevice, tableDevice);
    initializeExtrinsicInformation(double * L_e);

    for (int iter = 0; iter<niter; iter++) {
        
        deInterLeave<<<1,L_TOTAL>>>(L_eDevice, L_aDevice, tableDevice);

        gammaAlpha<<<1,L_TOTAL>>>(msgDevice , parity0Device,  L_aDevice,  gammaAlphaDevice);
        gammaBeta<<<1,L_TOTAL>>>(msgDevice , parity0Device,  L_aDevice,  gammaBetaDevice);
        cudaMemcpy(gammaAlphaHost, gammaAlphaDevice, sizeof(double)*L_TOTAL*4*4, cudaMemcpyDeviceToHost);
        cudaMemcpy(gammaBetaHost, gammaAlphaDevice, sizeof(double)*L_TOTAL*4*4, cudaMemcpyDeviceToHost);

        Alpha(AlphaHost, gammaAlphaHost);
        Beta(BetaHost, gammaBetaHost);
        cudaMemcpy(AlphaDevice, AlphaHost, sizeof(double)*L_TOTAL*4, cudaMemcpyHostToDevice);
        cudaMemcpy(BetaDevice, BetaHost, sizeof(double)*L_TOTAL*4, cudaMemcpyHostToDevice);
        normalizationAlphaAndBeta<<<1,L_TOTAL>>>(AlphaDevice, BetaDevice);

        LLRS(msgDevice, parity0Device, L_aDevice, AlphaDevice, BetaDevice, L_allDevice);

        extrinsicInformation<<<1, L_TOTAL>>>(L_allDevice, msgDevice, L_aDevice, L_eDevice);
        interLeave<<<1, L_TOTAL>>>(L_eDevice, L_aDevice, tableDevice);

        gammaAlpha<<<1,L_TOTAL>>>(imsgDevice , parity1Device,  L_aDevice,  gammaAlphaDevice);
        gammaBeta<<<1,L_TOTAL>>>(imsgDevice , parity1Device,  L_aDevice,  gammaBetaDevice);
        cudaMemcpy(gammaAlphaHost, gammaAlphaDevice, sizeof(double)*L_TOTAL*4*4, cudaMemcpyDeviceToHost);
        cudaMemcpy(gammaBetaHost, gammaAlphaDevice, sizeof(double)*L_TOTAL*4*4, cudaMemcpyDeviceToHost);

        Alpha(AlphaHost, gammaAlphaHost);
        Beta(BetaHost, gammaBetaHost);
        cudaMemcpy(AlphaDevice, AlphaHost, sizeof(double)*L_TOTAL*4, cudaMemcpyHostToDevice);
        cudaMemcpy(BetaDevice, BetaHost, sizeof(double)*L_TOTAL*4, cudaMemcpyHostToDevice);
        normalizationAlphaAndBeta<<<1,L_TOTAL>>>(AlphaDevice, BetaDevice);

        LLRS(imsgDevice, parity1Device, L_aDevice, AlphaDevice, BetaDevice, L_allDevice);

        extrinsicInformation<<<1, L_TOTAL>>>(L_allDevice, imsgDevice, L_aDevice, L_eDevice);
    }


    double L_all[L_TOTAL];
    undigned int msghat[L_TOTAL];

    // estimate information bits
    exestimateInformationBits(L_allDevice, mhatDevice, tableDevice); 
    cudaMemcpy(mhat, mhatDevice, sizeof(BYTE)*L_TOTAL, cudaMemcpyDeviceToHost);

    // count errors
    UINT bits_err = 0;
    for (i=0;i<L_TOTAL-M;i++) {
        if (mhat[i]!=m[i]) {
            bits_err++;
        }
    }
    cout<<"bits_err: "<<bits_err<<endl;

	delete m;
	delete x;
	delete y;
	delete mhat;

	return 0;
}
