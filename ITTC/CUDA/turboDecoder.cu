#include "helper_cuda.h"
#include "main.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
using namespace std;

#define M	3	// register length,=tail length
#define L_TOTAL 6144// if u want to use block interleave,L_TOTAL must = x^2
#define L_TOTAL_NUM 6147 
#define MAXITER 15
#define	FRAME_NUM 10000
//#define AlphaBetaBLOCK_NUM 8
//#define AlphaBetaTHREAD_NUM 8

//#define THREAD_NUM 8
#define BLOCK_NUM 12
#define L_BLOCK (L_TOTAL/BLOCK_NUM/4)
dim3 gridSize(2, BLOCK_NUM);
dim3 blockSize(4, 8);

#define LEAVER_BLOCK 8
#define LEAVER_THREAD 768
#define NSTATE	8	// = M^2
#define L_ALL (3*L_TOTAL+4*M)	// coded frame length

/*==================================================*/
/* ������� */
#define TERMINATED          1               /* 0����β��1��β */

#define TYPE_DECODER		1				/* ����������:1-LogMAP����
														  2-MAX-LogMAP����
														  3-SOVA����
                                                          4-const_LogMAP */
//#define N_ITERATION			8				/* ����������� */
#define MAX_FRAME_LENGTH	10000			/* ���֡�� */
/*==================================================*/
/* ��������� */
#define COLUMN_OF_G		4					/* ���������� */
#define G_ROW_1			13                  /* ������ͷ */					 
#define G_ROW_2			15					/* �����ͷ */

/*==================================================*/
/*==================================================*/

/* �������ݽṹ */

/* ������ṹ */
typedef struct
{
	int N_num_row;							/* ���������� */
	int K_num_col;							/* ���������� */
	int *g_matrix;							/* ��������ַ */
} TURBO_G;

/* Trellis�ṹ */
/*
mx_nextout(mx_lastout)Ϊ��(ǰ)�����������ַ:
����:״̬��, ����:4
ÿ�е�һ�к͵�����Ϊ��(ǰ)�������(1��-1),�ڶ��к͵�����Ϊ��֮��Ӧ�����(1��-1).

mx_nextstat(mx_laststat)Ϊ��(ǰ)��״̬������ַ:
����:״̬��, ����:2
���б�ʾ����Ϊ1(0)ʱ��Ӧ�ĺ�(ǰ)��״̬.
*/
typedef struct
{
	int *mx_nextout;		/* ����������� */						
	int *mx_nextstat;		/* ����״̬���� */
	int *mx_lastout;		/* ǰ��������� */
	int *mx_laststat;		/* ǰ��״̬���� */
	
} TURBO_TRELLIS;

/*==================================================*/
/* ����� */
#ifndef INFTY
#define INFTY 1E20
#endif
/*==================================================*/

/* ȫ�ֱ��� */

int *index_randomintlvr;		/* ҵ����Ϣ�����֯���±� */

int *index_randomintlvr2;

TURBO_G turbo_g;				/* ������ */

TURBO_TRELLIS turbo_trellis;	/* Tellis�ṹ */	

double rate_coding;				/* �������� */

void gen_qpp_index(int length, int *index);
/*==================================================*/
/* ����������صĲ��� */
int M_num_reg = COLUMN_OF_G-1;		/* �Ĵ����� */
int n_states = 8;						/* ״̬��:2��M_num_reg���� */
/*==================================================*/

int length_after_code;
int source_length;
int MODULATION;
int SYMBOL_NUM;
int f1,f2;
double _bpsk_map_i[2]=
{
	-1,1
};

double _bpsk_map_q[2]=
{
	0,0
};

__device__ float maxL(float x, float y) {
	return x>y?x:y;
}

__device__ float maxArray(float* arr, UINT length) {
	float temp = arr[0];
	for (int i = 1; i< length; i++) {
		if (arr[i] > temp)
			temp = arr[i];
	}
	return temp;
}
/*---------------------------------------------------------------
����:
	double E_algorithm(double x, double y)
����:
	Log-MAP�е�E�㷨:log(exp(x) + exp(y)) = max(x,y)+f(-|y-x|).
����:
	�������:
		x,y - ������.
	�������:
		��
����ֵ:
	�����
---------------------------------------------------------------*/
__device__ float E_algorithm(float x, float y)
{
	const double lookup_index_Log_MAP[16] = {0.0, 0.08824, 0.19587, 0.31026, 0.43275, 0.56508,
								0.70963, 0.86972, 1.0502, 1.2587, 1.5078, 1.8212,
								2.2522, 2.9706, 3.6764, 4.3758};
	const double lookup_table_Log_MAP[16] = {0.69315, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35,
								0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.025, 0.0125};
	float temp = (y-x)>0? (y-x):(x-y);
	int i;

	if (temp>=4.3758)
	{
		temp = 0;
	}
	else
	{
		/* �����±� */
		for (i=0; i<16 && temp>=lookup_index_Log_MAP[i]; i++)
		{
			;
		}
		/* ����f(-|y-x|) */
		temp = (float)lookup_table_Log_MAP[i-1];
	}
	
	/* ����max(x,y)+f(-|y-x|) */
	return ( (x>y?x:y) + temp );
}	

/*---------------------------------------------------------------
����:
	double E_algorithm_seq(double *data_seq, int length)
����:
	���е�E�㷨.
����:
	�������:
		data_seq - ������ַ.
		length - ���г���
	�������:
		��
����ֵ:
	�����
---------------------------------------------------------------*/

__device__ float E_algorithm_seq(float *data_seq, int length)
{
	int i;			/* ѭ������ */
	float temp;

	/* ÿ��������E�㷨,������һ������E�㷨 */
	temp = E_algorithm(*(data_seq+0), *(data_seq+1));
	for (i=2; i<length; i++)
	{
		temp = E_algorithm(temp, *(data_seq+i));
	}
	return temp;
}
//////////////////////////////////////////////////////////////////////
// LogMAP component decoder
// index true decoder1 false decoder2
//////////////////////////////////////////////////////////////////////
__global__ void logmap(float *msg, float* parity, float* L_a, float* L_all)
{
    const char NextOut[2][NSTATE] = // check bit based on current and input bit
    {	-1,-1,1,1,1,1,-1,-1,
        1,1,-1,-1,-1,-1,1,1
    };
    // NextState[bk][current state]
    const BYTE NextState[2][NSTATE] = // next state based on current and input bit
    {	0,4,5,1,2,6,7,3,
        4,0,1,5,6,2,3,7
    };
    // LastOut[bk][current state]
    //const char LastOut[2][NSTATE] =	// trellis last check bit
    //{	-1,1,1,-1,-1,1,1,-1,
    //    1,-1,-1,1,1,-1,-1,1
    //};
    // LastState[bk][current state]
    const BYTE LastState[2][NSTATE] =	// last state lead to current state by input bk
    {	0,3,4,7,1,2,5,6,
        1,2,5,6,0,3,4,7
    };

    //const unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    const unsigned int block = 4*(blockIdx.x*BLOCK_NUM + blockIdx.y) + threadIdx.x;
	//const unsigned int decoderIndex = unsigned int(block/BLOCK_NUM);
    const unsigned int threadX = threadIdx.x;
	//const unsigned int blockNum = (unsigned int)(threadInBlock/8);
	const unsigned int threadY = threadIdx.y;
	const unsigned int kIndex = blockIdx.x*(6144+3) + (block%(BLOCK_NUM*4))*L_BLOCK;

	float gamma0, gamma1;

	INT k;

	__shared__ float Alpha[L_BLOCK+4][4][8];
	__shared__ float Beta[2][4][8];

	__shared__ float tempSum0[4][8];
	__shared__ float tempSum1[4][8];

	// initialize Alpha & Beta
	if ((block == 0 || block == BLOCK_NUM*4)&& threadY != 0) {
			Alpha[0][threadX][threadY]=-INFTY;
	}
	else {
			Alpha[0][threadX][threadY] = 0;
	}

	// forward recursion,compute Alpha 
	for (k=1;k<L_BLOCK+4;k++)
	{
        gamma0=-msg[kIndex + k-1]+parity[kIndex + k-1]*NextOut[0][LastState[0][threadY]]
            -L_a[kIndex + k-1]/2;
        gamma1=msg[kIndex + k-1]+parity[kIndex + k-1]*NextOut[1][LastState[1][threadY]]
            +L_a[kIndex + k-1]/2;

		Alpha[k][threadX][threadY] = 
			maxL(gamma0 + Alpha[k-1][threadX][LastState[0][threadY]], 
				gamma1 + Alpha[k-1][threadX][LastState[1][threadY]]);
	}

	// backward recursion,compute Beta
    if ((block == BLOCK_NUM*4 - 1 || block == BLOCK_NUM*8 -1) && threadY != 0){
        Beta[1][threadX][threadY] = -INFTY;
    }
    else
        Beta[1][threadX][threadY] = 0;

    if (block == BLOCK_NUM*4 - 1 || block == BLOCK_NUM*8 -1){
		for (k=L_BLOCK+2;k>L_BLOCK-1;k--){
			gamma0=-msg[kIndex + k]+parity[kIndex + k]*NextOut[0][threadY]
				-L_a[kIndex + k]/2;
			gamma1=msg[kIndex + k]+parity[kIndex + k]*NextOut[1][threadY]
				+L_a[kIndex + k]/2;

			Beta[0][threadX][threadY] = 
				maxL(gamma0 + Beta[1][threadX][NextState[0][threadY]], 
					gamma1 + Beta[1][threadX][NextState[1][threadY]]);

			__syncthreads();
			gamma0=-msg[kIndex + k]+parity[kIndex + k]*NextOut[0][LastState[0][threadY]]
				-L_a[kIndex + k]/2;
			gamma1=msg[kIndex + k]+parity[kIndex + k]*NextOut[1][LastState[1][threadY]]
				+L_a[kIndex + k]/2;

        	tempSum0[threadX][threadY] = gamma0+Alpha[k][threadX][LastState[0][threadY]]+Beta[1][threadX][threadY];

        	tempSum1[threadX][threadY] = gamma1+Alpha[k][threadX][LastState[1][threadY]]+Beta[1][threadX][threadY];

			Beta[1][threadX][threadY]=Beta[0][threadX][threadY];
        	__syncthreads();

			if (threadY == 0) {
				L_all[kIndex + k]= maxArray(*(tempSum1+threadX), 8) - maxArray(*(tempSum0+threadX), 8); 
			}
		}
    } 

	for (k=L_BLOCK-1;k>=0;k--)
	{
		gamma0=-msg[kIndex + k]+parity[kIndex + k]*NextOut[0][threadY]
			-L_a[kIndex + k]/2;
		gamma1=msg[kIndex + k]+parity[kIndex + k]*NextOut[1][threadY]
			+L_a[kIndex + k]/2;

		Beta[0][threadX][threadY] = 
			maxL(gamma0 + Beta[1][threadX][NextState[0][threadY]], 
				gamma1 + Beta[1][threadX][NextState[1][threadY]]);

		__syncthreads();

		gamma0=-msg[kIndex + k]+parity[kIndex + k]*NextOut[0][LastState[0][threadY]]
			-L_a[kIndex + k]/2;
		gamma1=msg[kIndex + k]+parity[kIndex + k]*NextOut[1][LastState[1][threadY]]
			+L_a[kIndex + k]/2;

        tempSum0[threadX][threadY] = gamma0+Alpha[k][threadX][LastState[0][threadY]]+Beta[1][threadX][threadY];
        tempSum1[threadX][threadY] = gamma1+Alpha[k][threadX][LastState[1][threadY]]+Beta[1][threadX][threadY];

		Beta[1][threadX][threadY]=Beta[0][threadX][threadY];
        __syncthreads();

        if (threadY == 0) {
            L_all[kIndex + k]= maxArray(*(tempSum1+threadX), 8) - maxArray(*(tempSum0+threadX), 8); 
        }
	}
}

//__global__ void interLeave(float * src, float * des){
//    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
//    des[tid] = src[(((263 + tid*480)%6144)*tid)%6144];
//    des[(((263 + tid*480)%6144)*tid)%6144 + 6144] = src[tid+6144];
//}

__global__ void deInterLeave(float * src, float * des){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    //des[interLeaveTable[tid]] = src[tid];
    des[(((263 + tid*480)%6144)*tid)%6144] = src[tid];
    des[tid + 6144+3] = src[(((263 + tid*480)%6144)*tid)%6144 + 6144 +3];

	if (tid == 0){
	for (int i = 0; i < M; i++){
	    des[6144+i] = 0;
	    des[6144+3+6144+i] = 0;
	    }
	}
}

__global__ void extrinsicInformation(float * L_all, float * msg, float * L_a, float * L_e) {
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    L_e[tid] = L_all[tid + 6144+3] - 2*msg[tid + 6144+3] - L_a[tid + 6144+3];
    L_e[tid + 6144+3] = L_all[tid] - 2*msg[tid] - L_a[tid];

    if (tid == 0){
		for (int i = 0; i < 3; i++){
			L_e[6144+i] = L_all[6144+3+6144+i] - 2*msg[6144+3+6144+i] - L_a[6144+3+6144+i];
			L_e[6144 + 6144+3+i] = L_all[6144+i] - 2*msg[6144+i] - L_a[6144+i];
		}
    }
}

__global__ void demultiplex(float * stream, float * msg, float * parity) {
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    //if (puncture){// punctured rate=1/2
    //    msg[tid]=stream[2*tid];
    //    parity[tid%2][tid]=stream[tid*2+1];
    //}
    //else {// unpunctured rate=1/3
    //    msg[tid]=stream[3*tid];
    //    parity0[tid]=stream[3*tid+1];
    //    parity1[tid]=stream[3*tid+2];
    //}
        msg[tid]=0.5*stream[3*tid];
		msg[6144 + 3 + tid] = 0.5*stream[3*((((263 + tid*480)%6144)*tid)%6144)];
        parity[tid]=0.5*stream[3*tid+1];
        parity[6144+3+tid]=0.5*stream[3*tid+2];

		if (tid == 0){
		    for (int i = 0; i<M; i++){
		        msg[6144+i] = 0.5*stream[3*6144+2*i];
		        parity[6144+i] = 0.5*stream[3*6144+2*i+1];
				msg[6144 + 3+6144 + i] = 0.5*stream[3*6144+2*3+2*i];
				parity[6144 + 3+6144 + i] = 0.5*stream[3*6144+2*3+2*i+1];
			}
		}
}

__global__ void initializeExtrinsicInformation(float * L_e) {
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    L_e[tid] = 0;
	L_e[6144+3+tid] = 0;

	if (tid == 0){
	    for (int i = 0; i < M; i++){
	        L_e[6144+i] = 0;
	        L_e[6144+3+6144+i] = 0;
		}
	}
}

__global__ void exestimateInformationBits(float * L_all, int * msghat) {
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if(L_all[tid + 6144+3]>0)
        //msghat[m_Inter_table[tid]]=1;
		msghat[(((263 + tid*480)%6144)*tid)%6144] = 1;
    else
        //msghat[m_Inter_table[tid]]=0;
		msghat[(((263 + tid*480)%6144)*tid)%6144] = 0;
}

void countErrors(int *m, int * mhat, UINT * bitsError, UINT * frameError, UINT iter) {

	bool f_err = false;
	for (int i=0; i<(L_TOTAL);i++) {
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
	MODULATION = 1;			//���ƽ�����1,2,3,4,6
	source_length = 6144;			//������Դ����
	length_after_code = 6144*3+12;//SYMBOL_NUM*MODULATION;	//������׺�ĳ���	
	//SYMBOL_NUM=length_after_code/MODULATION;//  atoi(argv[4]);//672
	SYMBOL_NUM=length_after_code;

	f1 =263;
	f2 =480; 					//���볤�йصĽ�֯������

	double EbN0start = 0;//0;		//������ʼ�����
	double EbN0end = 1;//4;		//��������ֹ�����
	double EbN0step =0.1;		//��������Ȳ���
		
	double rate = (double)source_length/(double)(SYMBOL_NUM);

	int *source = NULL;
	int *mhat = NULL;
	clock_t start, finish;
	long double duration[FRAME_NUM];

	UINT bits_all,bits_err[MAXITER],frame_err[MAXITER];

	float Ber,Fer;

	int *coded_source = NULL;

	double *modulated_source_i, *modulated_source_q;

	double *after_channel_i,*after_channel_q;
	float *flow_for_decode = NULL;
	//int *flow_decoded = NULL;

	double EbN0dB,sigma;
	int nf, i1=0;

	TurboCodingInit();
			
	if ((source=(int *)malloc(source_length*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of source \n");
	  exit(1);  
	}
	mhat=(int *)malloc(source_length*sizeof(int));

	if ((coded_source=(int *)malloc((3*source_length+4*M_num_reg)*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of coded_source \n");
	  exit(1);  
	}

	if ((modulated_source_i=(double *)malloc((int)(SYMBOL_NUM)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of modulated_source_i \n");
	  exit(1);  
	}
	if ((modulated_source_q=(double *)malloc((int)(SYMBOL_NUM)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of modulated_source_q \n");
	  exit(1);  
	}
	if ((after_channel_i=(double *)malloc((int)(SYMBOL_NUM)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of after_channel_i \n");
	  exit(1);  
	}
	if ((after_channel_q=(double *)malloc((int)(SYMBOL_NUM)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of after_channel_q \n");
	  exit(1);  
	}

	if ((flow_for_decode=(float *)malloc((3*source_length+4*M_num_reg)*sizeof(float)))==NULL)
	{
	  printf("\n fail to allocate memory of flow_for_decode \n");
	  exit(1);  
	}
	//if ((flow_decoded=(int *)malloc(N_ITERATION*source_length*sizeof(int)))==NULL)
	//{
	//  printf("\n fail to allocate memory of flow_decoded\n");
	//  exit(1);  
	//}
	
	srand((unsigned)time(NULL));

	findCudaDevice(argc, (const char **)argv);

	float * yDevice;
	float * msgDevice;
	int * mhatDevice;
	float * parityDevice;
	float * L_eDevice;
	float * L_aDevice;
	float * L_allDevice;

    cudaMalloc((void **)&yDevice, L_ALL*sizeof(float));
    cudaMalloc((void **)&msgDevice, L_TOTAL_NUM*2*sizeof(float));
    cudaMalloc((void **)&mhatDevice, L_TOTAL*sizeof(int));
    cudaMalloc((void **)&parityDevice, L_TOTAL_NUM*2*sizeof(float));
    cudaMalloc((void **)&L_eDevice, L_TOTAL_NUM*2*sizeof(float));
    cudaMalloc((void **)&L_aDevice, L_TOTAL_NUM*2*sizeof(float));
    cudaMalloc((void **)&L_allDevice, L_TOTAL_NUM*2*sizeof(float));

	for (EbN0dB=EbN0start; EbN0dB<=EbN0end; EbN0dB+=EbN0step)
	{
		sigma = pow(10,-EbN0dB/20)*sqrt(0.5/(rate*MODULATION));

		for (int i =0; i<MAXITER;i++) {
			bits_err[i]=0;
			frame_err[i]=0;
		}
		bits_all = 0;

		for (nf=0; nf<FRAME_NUM; nf++,bits_all += L_TOTAL)
		{
			for(i1=0; i1<source_length; i1+=1)
			{
				*(source+i1)=rand()%2;
			}

            TurboEnCoding(source, coded_source, source_length);

/*******************************************************************************/
			module(coded_source,modulated_source_i,modulated_source_q,SYMBOL_NUM*MODULATION,MODULATION);

			AWGN(modulated_source_i, after_channel_i, sigma, SYMBOL_NUM);
			AWGN(modulated_source_q, after_channel_q, sigma, SYMBOL_NUM);

			demodule(after_channel_i, after_channel_q, SYMBOL_NUM,flow_for_decode,1/(2*pow(sigma,2)),MODULATION);
/*******************************************************************************/

			cudaMemcpy(yDevice,flow_for_decode,sizeof(float)*length_after_code, cudaMemcpyHostToDevice);
			start = clock();

			demultiplex<<<LEAVER_BLOCK,LEAVER_THREAD>>>(yDevice, msgDevice, parityDevice); 
			initializeExtrinsicInformation<<<LEAVER_BLOCK,LEAVER_THREAD>>>(L_eDevice);

			for (int iter = 0; iter<MAXITER; iter++) {
				
				deInterLeave<<<LEAVER_BLOCK,LEAVER_THREAD>>>(L_eDevice, L_aDevice);

                logmap<<<gridSize, blockSize>>>(msgDevice, parityDevice, L_aDevice, L_allDevice);

				extrinsicInformation<<<LEAVER_BLOCK,LEAVER_THREAD>>>(L_allDevice, msgDevice, L_aDevice, L_eDevice);

				//interLeave<<<LEAVER_BLOCK,LEAVER_THREAD>>>(L_eDevice, L_aDevice);

                //logmap<<<BLOCK_NUM, THREAD_NUM>>>(imsgDevice, parity1Device, L_aDevice, L_allDevice, false);

				//extrinsicInformation<<<LEAVER_BLOCK,LEAVER_THREAD>>>(L_allDevice, imsgDevice, L_aDevice, L_eDevice);

				exestimateInformationBits<<<LEAVER_BLOCK,LEAVER_THREAD>>>(L_allDevice, mhatDevice); 

				cudaMemcpy(mhat, mhatDevice, sizeof(int)*L_TOTAL, cudaMemcpyDeviceToHost);
				countErrors(source, mhat, bits_err, frame_err, iter);

			}
			finish = clock();
			duration[nf] = (long double)(finish - start);

			//cudaMemcpy(mhat, mhatDevice, sizeof(int)*L_TOTAL, cudaMemcpyDeviceToHost);
			//countErrors(source, mhat, bits_err, frame_err, MAXITER-1);

		}
		printf("-------------------------\n");
		printf("Eb/No=%fdB:\n",EbN0dB);
		printf("-------------------------\n");

		for (int i=0;i<MAXITER;i++) {
			Ber=(float)bits_err[i]/(float)bits_all;
			Fer=(float)frame_err[i]/(float)FRAME_NUM;
			printf("Iteration:%d\n",i+1);
			printf("---Ber=%f\n---Fer=%f\n",Ber,Fer);
		}
		long double durationSum = 0.0;
		for (int i = 0; i < FRAME_NUM; i++) {
			durationSum += duration[i];
		}
		durationSum /= CLOCKS_PER_SEC;

		long double throughput = FRAME_NUM*(6144+3) / durationSum / 1000000;
		cout<<"throughput: "<<throughput<<"Mbps"<<endl;
			
	}
	
/*-----------------------------------------------------------------*/	
/*-----------------------------------------------------------------*/	
	TurboCodingRelease();
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/	

	free(source);
	free(coded_source);
	free(modulated_source_i);
	free(modulated_source_q);
	free(after_channel_i);
	free(after_channel_q);
	free(flow_for_decode);
	//free(flow_decoded);

	cudaFree(yDevice);
	cudaFree(msgDevice);
	cudaFree(mhatDevice);
	cudaFree(parityDevice);
	cudaFree(L_eDevice);
	cudaFree(L_aDevice);
	cudaFree(L_allDevice);
} //main


void randominterleaver_long(long *data_unintlvr, long *interleaverddata, int *index_randomintlvr, int length)
{
	int i;
	int *index_random = index_randomintlvr;

	for (i=0; i<length; i++)
	{
		*(interleaverddata+i) = *(data_unintlvr+ (*(index_random+i)));
	}
}

void random_deinterlvr_long(long *data_unintlvr, long *interleaverddata, int *index_randomintlvr, int length)
{
	int i;
	int *index_random = index_randomintlvr;

	for (i=0; i<length; i++)
	{
		*(data_unintlvr+(*(index_random+i))) = *(interleaverddata+i);
	}
}


void randominterleaver_int(int *data_unintlvr, int *interleaverddata, int *index_randomintlvr, int length)
{
	int i;
	int *index_random = index_randomintlvr;

	for (i=0; i<length; i++)
	{
		*(interleaverddata+i) = *(data_unintlvr+ (*(index_random+i)));
	}
}

void random_deinterlvr_int(int *data_unintlvr, int *interleaverddata, int *index_randomintlvr, int length)
{
	int i;
	int *index_random = index_randomintlvr;

	for (i=0; i<length; i++)
	{
		*(data_unintlvr+(*(index_random+i))) = *(interleaverddata+i);
	}
}

void randominterleaver_double(double *data_unintlvr, double *interleaverddata, int *index_randomintlvr, int length)
{
	int i;
	int *index_random = index_randomintlvr;

	for (i=0; i<length; i++)
	{
		*(interleaverddata+i) = *(data_unintlvr+ (*(index_random+i)));
	}
}

void random_deinterlvr_double(double *data_unintlvr, double *interleaverddata, int *index_randomintlvr, int length)
{
	int i;
	int *index_random = index_randomintlvr;

	for (i=0; i<length; i++)
	{
		*(data_unintlvr+(*(index_random+i))) = *(interleaverddata+i);
	}
}

/*---------------------------------------------------------------
����:
	int gen_g_matrix(int k_column, int g_row1, int g_row2, int *mx_g_turbo)
����:
	�õ�������.
����:
	�������:
		k_column - ����������.
		g_row1 - �������һ��.
		g_row2 - ������ڶ���.
	�������:
		mx_g_turbo - ��������ַ.
����ֵ:
	1 - �ɹ��صõ�������.
	0 - ����������ʧ��.
---------------------------------------------------------------*/
int gen_g_matrix(int k_column, int g_row1, int g_row2, int *mx_g_turbo)
{
	int i, position;		/* ѭ������ */
	int high_num, low_num;
	
	/* ��һ�� */
	high_num = g_row1;		
	position = 1;			/* �ڼ���8������ */
	while (high_num>0)
	{
		low_num = high_num%10;	/* �õ����8������ */
		if (low_num>7)			/* �ж��Ƿ�Ϊ8������ */
		{
			return 0;
		}
		high_num = high_num/10;		/* �������ಿ�� */

		/* ��8������תΪ�����Ʋ����� */
		for (i=k_column-(position-1)*3-1; i>=0 && i>=k_column-position*3; i--)
		{
			*(mx_g_turbo+i) = low_num%2;
			low_num = low_num/2;
		}
		position++;		/* ��һ��λ�� */
		if (i<0)
		{
			break;
		}
	}

	/* �ڶ��� */
	high_num = g_row2;
	position = 1;			/* �ڼ���8������ */
	while (high_num>0)
	{
		low_num = high_num%10;		/* �õ����8������ */
		if (low_num>7)				/* �ж��Ƿ�Ϊ8������ */
		{
			return 0;
		}
		high_num = high_num/10;		/* �������ಿ�� */

		/* ��8������תΪ�����Ʋ����� */
		for (i=k_column-(position-1)*3-1; i>=0 && i>=k_column-position*3; i--)
		{
			*(mx_g_turbo+k_column+i) = low_num%2;
			low_num = low_num/2;
		}
		position++;					/* ��һ��λ�� */
		if (i<0)
		{
			break;
		}
	}
	return 1;
}

/*---------------------------------------------------------------
����:
	void int2bin(int intstat, int *tempstat, int length)
����:
	ʮ������תΪ����������.
����:
	�������:
		intstat - ʮ������.
		length - Ҫ�õ��Ķ��������г��ȣ�
	�������:
		bin_stat - ������������ַ.
����ֵ:

	��
---------------------------------------------------------------*/
void int2bin(int intstat, int *bin_stat, int length)
{
	int i, temp;

	temp = intstat;

	/* ����2������ */
	for (i=length-1; i>=0; i--)
	{
		*(bin_stat+i) = temp%2;
		temp = temp/2;
	}
}

/*---------------------------------------------------------------
����:
	int bin2int(int *binseq, int length)
����:
	����������תΪʮ������.
����:
	�������:
		binseq - ������������ַ.
		length - ���������г��ȣ�
	�������:
		��
����ֵ:
	�õ���ʮ��������
---------------------------------------------------------------*/
int bin2int(int *binseq, int length)
{
	int i, j, temp;
	int sum = 0;

	for (i=0; i<length; i++)
	{
		temp = 1;

		/* �����λȨֵ */
		for (j=1; j<=i; j++)
		{
			temp = temp * 2;
		}
		sum = sum + temp * (*(binseq+length-1-i));
	}

	return sum;
}

/*---------------------------------------------------------------
����:
	int encode_bit(int inbit, int *stat)
����:
	���ر�����.
����:
	�������:
		inbit -�������.
		stat - ��ǰ�Ĵ���״̬��ַ.
	�������:
		stat - �����Ĵ���״̬��ַ.
����ֵ:
	������أ�
---------------------------------------------------------------*/
int encode_bit(int inbit, int *stat)
{
	int j;			/* ѭ������ */
	int output;		/* ������� */

	/* ����������� */
	output = (*(turbo_g.g_matrix+turbo_g.K_num_col+0)) * inbit;

	for (j=1; j<turbo_g.K_num_col; j++)
	{
		output = (output + (*(turbo_g.g_matrix+turbo_g.K_num_col+j)) * (*(stat+j-1)))%2;
	}

	/* �޸�״̬���� */
	for (j=turbo_g.K_num_col-2; j>0; j--)
	{
		*(stat+j)=*(stat+j-1);
	}

	*(stat+0) = inbit;

	return output;
}

/*---------------------------------------------------------------
����:
	void gen_trellis()
����:
	����Trellis.
����:
	��
����ֵ:
	��
---------------------------------------------------------------*/
void gen_trellis()
{
	int i, j, k;						/* ѭ������ */
	int dk_turbo, ak_turbo, outbit;		/* �������ڲ����غ�������� */

	int *tempstat;						/* ״̬���� */

	if ((tempstat=(int *)malloc(sizeof(int)*M_num_reg))==NULL)
	{
	  printf("\n fail to allocate memory of tempstat \n");
	  exit(1);  
	}
	/* ���ɺ�������ͺ���״̬���� */
	for (i=0; i<n_states; i++)			/* ״̬ѭ�� */
	{
		for (j=0; j<2; j++)				/* ����Ϊ0,1 */
		{
			int2bin(i, tempstat, M_num_reg);	/* ��״̬תΪ���������� */

			/* dk:ԭʼ������� */
			dk_turbo = j;

			/* ����ak:���ӷ������������أ��������ɾ����һ�У�13�� */
			ak_turbo = (*(turbo_g.g_matrix+0)) * dk_turbo;
			for (k=1; k<turbo_g.K_num_col; k++)
			{
				ak_turbo = ak_turbo + (*(turbo_g.g_matrix+k)) * (*(tempstat+k-1));
			}

			ak_turbo = ak_turbo % 2;

			/* �����������,�޸�״̬���У��������ɾ���ڶ��У�15�� */
			outbit = encode_bit(ak_turbo, tempstat);

			/* д���������ͺ���״̬���� */
			*(turbo_trellis.mx_nextout+i*4+2*j)=2*dk_turbo-1;
			*(turbo_trellis.mx_nextout+i*4+2*j+1)=2*outbit-1;
			*(turbo_trellis.mx_nextstat+i*2+j)=bin2int(tempstat, M_num_reg);
		}			/* ����ѭ������ */

	}	/* ״̬ѭ������ */

	/* ����ǰ�������ǰ��״̬���� */
	for (i=0; i<n_states; i++)	/* ״̬ѭ�� */
	{
		for (j=0; j<2; j++)/* ����Ϊ0,1 */
		{
			*(turbo_trellis.mx_laststat+(*(turbo_trellis.mx_nextstat+i*2+j))*2+j) = i;
			*(turbo_trellis.mx_lastout+(*(turbo_trellis.mx_nextstat+i*2+j))*4+2*j)
				= *(turbo_trellis.mx_nextout+i*4+2*j);
			*(turbo_trellis.mx_lastout+(*(turbo_trellis.mx_nextstat+i*2+j))*4+2*j+1)
				= *(turbo_trellis.mx_nextout+i*4+2*j+1);
		}	/* ����ѭ������ */
	}	/* ״̬ѭ������ */

	free(tempstat);
}

/*---------------------------------------------------------------
����:
	void TurboCodingInit()
����:
	Turbo���ʼ������.
����:
	��
����ֵ:
	��
---------------------------------------------------------------*/
void TurboCodingInit()
{
	/* ��ʼ�������� */
	turbo_g.N_num_row = 3;				/* ���� */
	turbo_g.K_num_col = COLUMN_OF_G;	/* ���� */

	SYMBOL_NUM = (int)((double)length_after_code/(MODULATION));

	if ((turbo_g.g_matrix=(int *)malloc(turbo_g.N_num_row*turbo_g.K_num_col*sizeof(int)))==NULL)
	{
		printf("\n fail to allocate memory of turbo_g\n");
		exit(1);  
	}

	/* �õ������� */
	if (!gen_g_matrix(turbo_g.K_num_col, G_ROW_1, G_ROW_2, turbo_g.g_matrix))
	{
		printf("error number of G\n");
		exit(1);
	}

	/* ����Trellis */
	if ((turbo_trellis.mx_lastout=(int *)malloc(sizeof(int)*n_states*4))==NULL)
	{
	  printf("\n fail to allocate memory of turbo_trellis.mx_lastout \n");
	  exit(1);  
	}
	if ((turbo_trellis.mx_laststat=(int *)malloc(sizeof(int)*n_states*2))==NULL)
	{
	  printf("\n fail to allocate memory of turbo_trellis.mx_laststat \n");
	  exit(1);  
	}
	if ((turbo_trellis.mx_nextout=(int *)malloc(sizeof(int)*n_states*4))==NULL)
	{
	  printf("\n fail to allocate memory of turbo_trellis.mx_nextout \n");
	  exit(1);  
	}
	if ((turbo_trellis.mx_nextstat=(int *)malloc(sizeof(int)*n_states*2))==NULL)
	{
	  printf("\n fail to allocate memory of turbo_trellis.mx_nextstat \n");
	  exit(1);  
	}

	gen_trellis();

	/* Ϊ�����֯�������ڴ� */
		if ((index_randomintlvr=(int *)malloc(MAX_FRAME_LENGTH*sizeof(int)))==NULL)
		{
			printf("\n fail to allocate memory of index_randomintlvr \n");
			exit(1);  
		}

		if ((index_randomintlvr2=(int *)malloc(2*MAX_FRAME_LENGTH*sizeof(int)))==NULL)
		{
			printf("\n fail to allocate memory of index_randomintlvr \n");
			exit(1);  
		}

	/*����QPP��֯���±�*/
	gen_qpp_index(source_length, index_randomintlvr);

	/*int i1,*num;

	if ((num=(int *)malloc(source_length*sizeof(int)))==NULL)
	{
		printf("\n fail to allocate memory of num \n");
		exit(1);  
	}

	for (i1=0;i1<source_length;i1++)
	{
		num[i1]=0;	
	}
	for (i1=0;i1<source_length;i1++)
	{
		num[index_randomintlvr[i1]]++;
			
	}

	for (i1=0;i1<source_length;i1++)
	{
	printf("%d ",num[i1]);
			
	}
    free(num);*/
}

/*---------------------------------------------------------------
����:
	void rsc_encode(int *source, int *rsc, int terminated, int len_info)
����:
	RSC������.
����:
	�������:
		source -Դ����������ַ.
		len_info - Դ�������г��ȣ�
		terminated - �Ƿ��β: 1-��β, 0-����β.
	�������:
		RSC - ���������������ַ.
����ֵ:
	�ޣ�
---------------------------------------------------------------*/
void rsc_encode(int *source, int *rsc, int terminated, int len_info)
{
	int i, j;			/* ѭ������ */

	int *state;			/* ״̬���� */
	int dk_turbo, ak_turbo, outbit;		/* �������ڵ�dk,ak�ͱ���������� */

	int len_total;						/* �ܳ��� */

	if ((state=(int *)malloc(M_num_reg*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of state \n");
	  exit(1);  
	}

	/* �����ܳ��� */
	len_total = len_info+M_num_reg;

	/* ��״̬Ϊ0 */
	for (i=0; i<M_num_reg; i++)
	{
		*(state+i) = 0;
	}

	for (i=0; i<len_total; i++)							/* ����ر��� */
	{
		if (!terminated || (terminated && i<len_info))	/* ����Ϣ���� */
		{
			dk_turbo = *(source+i);
		}
		else											/* ��βTrellis */
		{
			if (terminated && i>=len_info)
			{
				dk_turbo = 0;
				for (j=1; j<turbo_g.K_num_col; j++)
				{
					dk_turbo = dk_turbo + (*(turbo_g.g_matrix+j)) * (*(state+j-1));
				}
				dk_turbo = dk_turbo%2;
			}
		}

		/* ����ak */
		ak_turbo = *(turbo_g.g_matrix+0) * dk_turbo;
		for (j=1; j<turbo_g.K_num_col; j++)
		{
			ak_turbo = ak_turbo + (*(turbo_g.g_matrix+j))*(*(state+j-1));
		}

		ak_turbo = ak_turbo%2;

		/* ��ak���б��ر��� */
		outbit = encode_bit(ak_turbo, state);

		/* дdk��������� */
		*(rsc+2*i) = dk_turbo;
		*(rsc+2*i+1) = outbit;
	}				/* ����ر������ */

	free(state);
}

/*---------------------------------------------------------------
����:
	void encoderm_turbo(int *source, int *send_turbo, int len_info)

����:
	Turbo�������.
����:
	�������:
		source -Դ����������ַ.
		len_info - Դ�������г��ȣ�
		type_flow - ��Ϣ����: 1-ҵ����Ϣ, 0������Ϣ.
	�������:
		send_turbo - ������ƺ�����������ַ.
����ֵ:
	�ޣ�
---------------------------------------------------------------*/
void encoderm_turbo(int *source, int *send_turbo, int len_info)
{
	int i;									/* ѭ������ */
	int len_total = len_info + M_num_reg;	/* �ܳ��� */

	int *rsc1, *rsc2;		/* ����RSC����������� */
	
	int *input2;			/* RSC2������ */

	if ((rsc1=(int *)malloc(2*len_total*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of rsc1 \n");
	  exit(1);  
	}
	if ((rsc2=(int *)malloc(2*len_total*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of rsc2 \n");
	  exit(1);  
	}

	if ((input2=(int *)malloc(len_info*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of input2 \n");
	  exit(1);  
	}

	/* RSC1 */
	rsc_encode(source, rsc1, TERMINATED, len_info);

	/* ��֯��Դ���� */
	randominterleaver_int(source, input2, index_randomintlvr, len_info);

	/* RSC2 */
	rsc_encode(input2, rsc2, TERMINATED, len_info);

	/* ��Ϣλ���õ��� */
	for (i=0; i<len_info; i++)
	{
		*(send_turbo+3*i) = *(rsc1+2*i);
		*(send_turbo+3*i+1) = *(rsc1+2*i+1);
		*(send_turbo+3*i+2) = *(rsc2+2*i+1);
	}
	
	/* ��βλ���õ��� */
	for (i=0; i<2*M_num_reg; i++)
	{
		*(send_turbo+3*len_info+i) = *(rsc1+2*len_info+i);
		*(send_turbo+3*len_info+2*M_num_reg+i) = *(rsc2+2*len_info+i);
	}
	
	free(rsc1);
	free(rsc2);
	free(input2);
}


/*---------------------------------------------------------------
����:
	double random_turbo()
����:
	����0-1���ȷֲ������.
����:
	��
����ֵ:
	���ɵ�0-1���ȷֲ��������
---------------------------------------------------------------*/
double random_turbo()
{
	long z,k;
	static long s1 = 12345L;
	static long s2 = 1234546346L;

	k= s1 / 53668L;
	s1 = 40014L * (s1 - k*53668L) - k*12211L;
	if (s1<0)
	  s1 = s1 + 2147483563L;
	k = s2 / 52774;
	s2 = 40692L * (s2 - k*52774L) - k*3791L;
	if (s2<0)
        s2 = s2 + 2147483399L;
 	z=s1 - s2;
	if (z<1)
  	  z = z + 2147483562L;
	return (double) z / (double) 2147483563.0;
}

void gen_qpp_index(int length, int *index)
{
	int i;

	for (i=0; i<source_length; i++)
	{
		*(index+i) = (f1*i+(((f2*i)%source_length)*i)%source_length)%source_length; // ������ת����Ϊ�˷�ֹ��ֵ̫��������
	}
}

/*---------------------------------------------------------------
����:
	void gen_rand_index(int length, int type_flow)
����:
	���������֯�����±�.
����:
	�������:
		length - ��֯������.
		type_flow - ��Ϣ����: 1-ҵ����Ϣ, 0-������Ϣ.
	�������:
		��
����ֵ:
	��
---------------------------------------------------------------*/
void gen_rand_index(int length, int *index)
{
	int *index_random;			/* �������ɵ�0-1���ȷֲ���������� */
	double *tempindex;			/* ����ѡ�������֯��ַ */
	double tempmax;				/* ���ֵ */
	int selectedscr;			/* ѡ�е��±� */
	int i, j;					/* ѭ������ */

	if ((tempindex=(double *)malloc((length)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of tempindex \n");
	  exit(1);  
	}

	/* ����������ѡ�����֯�� */
	index_random = index;

	/* ���ɵ�0-1���ȷֲ���������� */
	for (i=0; i<length; i++)
	{
		*(tempindex+i) = random_turbo();
	}

	for (i=0; i<length; i++)	
	{
		/* �ҵ�tempindex�е����ֵ��Ӧ���±� */
		tempmax = 0.0;

		for (j=0; j<length; j++)
		{
			if (*(tempindex+j) >= tempmax)
			{
				tempmax = *(tempindex+j);
				selectedscr = j;
			}
		}

		/* ��֯����λ��Ϊ���±� */
		*(index_random+i) = selectedscr;

		/* tempindex��λ��0 */
		*(tempindex+selectedscr) = 0.0;
	}

	free(tempindex);
}

/*---------------------------------------------------------------
����:
	TurboEnCoding(int *source, int *coded_source,
						int *source_length)
����:
	Turbo��ҵ����Ϣ���뺯��.
����:
	�������: source - Դbit������ַ.
			  source_length - Դbit���г���.
	�������: coded_source - ������������ַ.
����ֵ:
	��
---------------------------------------------------------------*/
void TurboEnCoding(int *source, int *coded_source, int source_length)
{
	int i;							/* ѭ������ */

	int *temp_send = NULL;			
	int *send = NULL;				/* ����������ַ */

	int length_info = source_length;		/* ��Ϣλ���� */

	/* �����ڴ� */
	if ((send=(int *)malloc((3*length_info+4*M_num_reg)*sizeof(int)))==NULL)
	{
		printf("\n fail to allocate memory of send \n");
		exit(1);
	}

	/* ���������֯���±� */
//	gen_rand_index(length_info, index_randomintlvr);

	encoderm_turbo(source, send, length_info);	/* ���� */

	temp_send = send;

	/* д������� */
	for (i=0; i<(3*length_info+4*M_num_reg); i++)
	{
		*(coded_source+i) = *(temp_send+i);
	}

	free(send);
}

/*---------------------------------------------------------------
����:
	double get_max(double *data_seq, int length)
����:
	�õ����е����ֵ.
����:
	�������:
		data_seq - ������ַ.
		length - ���г���.
	�������:
		��
����ֵ:
	�����е����ֵ��
---------------------------------------------------------------*/
double get_max(double *data_seq, int length)
{
	int i;		/* ѭ������ */
	double temp;
	temp = *(data_seq+0);
	for (i=1; i<length; i++)
	{
		if (temp < *(data_seq+i))
		{
			temp = *(data_seq+i);
		}
	}

	return temp;
}







/*---------------------------------------------------------------
����:
	void decision(double *LLR_seq, int length, int *output)
����:
	�о�.
����:
	�������:
		LLR_seq - LLR������ַ.
		length - ���г���
	�������:
		output - �о������������ַ
����ֵ:
	�ޣ�
---------------------------------------------------------------*/
void decision(double *LLR_seq, int length, int *output)
{
	int i;

	for (i=0; i<length; i++)
	{
		/* С��0��0 */
		if (*(LLR_seq+i) < 0)
		{
			*(output+i) = 0;
		}
		/* ����0��1 */
		else
		{
			*(output+i) = 1;
		}
	}		
}


void dectobin(double *flow_for_change, int *flow_changed, int flow_len, int integer_len, int decimal_len)
{
    int i;/* ѭ������ */
	int total_len = integer_len + decimal_len;

		/* �����ڴ� */
	/*if ((flow_changed=(int *)malloc((flow_len)*sizeof(int)))==NULL)
	{
		printf("\n fail to allocate memory of receive_punc \n");
		exit(1);  
	}*/

	for (i=0; i<flow_len; i++)
	{
		*(flow_for_change+i) *= pow(2,decimal_len);
		*(flow_changed+i) = (int)(*(flow_for_change+i)+0.5);
		if (*(flow_changed+i)<-pow(2,total_len))
		{
			*(flow_changed+i)=(long)-pow(2,total_len);
		}
		else if (*(flow_changed)>(pow(2,total_len)-1))
		{
			*(flow_changed+i)=(long)pow(2,total_len+i)-1;
		}            //-__logf(1+__expf(L_a[block*L_BLOCK + k-1]));

		/*if (*(flow_changed+i)>=0)
		{
			*(flow_changed+i) = *(flow_changed+i);
		}
		else 
		{
			*(flow_changed+i) = *(flow_changed+i)+pow(2,total_len);
		}*/

	}
	
	/*ofstream flow_changed_f("flow_changed.txt");
	flow_changed_f<< "-----flow_changed�ļ���ʼ-----" << endl;
	for (i=0; i<flow_len; i++) flow_changed_f << flow_changed[i] << endl;

	flow_changed_f<< "-----flow_changed�ļ�����-----" << endl;*/
}


void TurboCodingRelease()
{
	/* �ͷ��������ڴ� */
	free(turbo_g.g_matrix);
	/* �ͷ�trellis�ṹ�ڴ� */
	free(turbo_trellis.mx_lastout);
	free(turbo_trellis.mx_laststat);
	free(turbo_trellis.mx_nextout);
	free(turbo_trellis.mx_nextstat);
	/* �ͷ������֯���ڴ� */
	free(index_randomintlvr);
	free(index_randomintlvr2);

}
/*---------------------------------------------------------------
����:
	void mgrns(double mean,double sigma,double seed,int n,double *a)
����:
	��������Ϊn�ĸ�˹�������.
����:
	�������:	mean - ��ֵ
				sigma - ��׼��
				seed - һ���������
	�������:	a - ����Ϊn�ĸ�˹�������.
����ֵ:
	��
---------------------------------------------------------------*/
void mgrns(double mean,double sigma,double seed,int n,double *a)
{ int i,k,m;
    double s,w,v,t;
    s=65536.0; w=2053.0; v=13849.0;
    for (k=0; k<=n-1; k++)
	{
		t=0.0;
		for (i=1; i<=12; i++)
        { 
			seed=seed*w+v; m=(int)(seed/s);
            seed=seed-m*s; t=t+(seed)/s;
        }/*�������ļ��޶���������Ӹ�˹�ֲ��������*/
        *(a+k)=mean+(double)(sigma*(t-6.0));
    }
    return;
}
/*---------------------------------------------------------------
����:
	void AWGN(int *send, double *r, double sigma, int totallength)
����:
	AWGN�ŵ�.
����:
	�������: send - int�ͷ���������ַ.
			  sigma - ������׼��.
			  totallength - ���г���.
	�������: r - ��AWGN�ŵ�����������е���ַ.
����ֵ:
	��
---------------------------------------------------------------*/
void AWGN(double *send, double *r, double sigma, int totallength)
{
	int i;
	double *noise = (double *)malloc(sizeof(double)*totallength);
	double seed =  (double)(3.0 - (double)((rand() & RAND_MAX)/(double)RAND_MAX)/10e6);
	mgrns(0,sigma,seed,totallength,noise);
	for(i=0; i<totallength; i++)
	{
		*(r+i) = (double)( *(send+i) + *(noise+i) );
		//*(r+i) = (double)( *(send+i) + 0 );
	}
	free(noise);
}

/************************************************************************/

double  calculate_sqr_dis(double pt1_r,double pt1_i,double pt2_r,double pt2_i)
{
	double dis_r,dis_i,sqr_dis;

	dis_r=pt1_r-pt2_r;
	
	dis_i=pt1_i-pt2_i;

	sqr_dis=dis_r*dis_r+dis_i*dis_i;

	return sqr_dis;


}

void _bpsk_module(int * a,double * outi,double * outq,int N)
{
	int i;

	for (i=0;i<N;i++)
	{
		*(outq+i)=0.0;
		
		if (*(a+i)==1) 
			*(outi+i)=1.0;
		else 
			*(outi+i)=-1.0;
	}

}


void module(int * a,double * outi,double * outq,int N, int modu_index)
{
	switch (modu_index)
	{
	    case 1: _bpsk_module(a,outi,outq,N);break;
	}

}


void _bpsk_demodule(double *symbol_i, double *symbol_q, int symbol_len,
					float *out,double Kf)
{
	int i,j;
	double min_sqr_dis1,min_sqr_dis2,sqr_dis1,sqr_dis2;

	for(i=0;i<symbol_len;i++)
	{
	/*the minimum square distance has nothing with the orders of calculation of sqr_dis1 and sqr_dis2*/

		min_sqr_dis1=0x7fffffffffff,min_sqr_dis2=0x7fffffffffff;

        for(j=0;j<2;j++)
		{
			if(j&1)
			{
			/*calculate Min{dj} where j in Sj*/
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_bpsk_map_i[j],(double)_bpsk_map_q[j]);
	
			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_bpsk_map_i[j],(double)_bpsk_map_q[j]);

			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}
		}
        out[i]= -Kf*(min_sqr_dis1-min_sqr_dis2);
	}
}

void demodule(double *symbol_i, double *symbol_q, int symbol_len,float* out,double Kf,int modu_index)
{
    switch(modu_index)
	{
	    case 1: _bpsk_demodule(symbol_i, symbol_q,symbol_len,out,Kf);break;
	}
}


