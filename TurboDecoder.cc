// TurboDecoder : Defines the entry point for the console application.
//
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

#define MAXITER 5
#define	FRAME_NUM 10

//void sova(double *msg, double *parity, double *L_a,double *L_all, bool index);
UINT m_Inter_table[L_TOTAL];





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

//////////////////////////////////////////////////////////////////////
// LogMAP component decoder
//////////////////////////////////////////////////////////////////////
void logmap(double *msg, double *parity, double *L_a,double *L_all, bool index)
{
	UINT nstate;
	UINT s1,s2;
	double (* Beta)[4];
	double (*Alpha)[4];
	double * max_branch;
	double gamma[4];
	double sum,sum0,sum1;
	double INIFINITY_ = 1E+10;	// approximate infinity value

	INT k;

	// alloc memory,
	Beta = new double[L_TOTAL+1][4];
	Alpha = new double[L_TOTAL+1][4];
	max_branch = new double[L_TOTAL+1];

	nstate=NSTATE;

	// initialize Alpha & Beta
	Alpha[0][0]=0;
	for (s1=1;s1<nstate;s1++)
		Alpha[0][s1]=-INIFINITY_;
	if (index){	// true -- terminated,false -- open
		Beta[L_TOTAL][0]=0;
		for (s2=1;s2<nstate;s2++)
			Beta[L_TOTAL][s2]=-INIFINITY_;
	}else{
		for (s2=0;s2<nstate;s2++)
			Beta[L_TOTAL][s2]=0;
	}

	// forward recursion,compute Alpha 
	for (k=1;k<=L_TOTAL;k++)
	{
		for (s2=0;s2<nstate;s2++)
		{
			sum=0;
			for (s1=0;s1<nstate;s1++)
				gamma[s1]=-INIFINITY_;
			gamma[LastState[0][s2]]=-msg[k-1]+parity[k-1]*LastOut[0][s2]
				-log(1+exp(L_a[k-1]));
			gamma[LastState[1][s2]]=msg[k-1]+parity[k-1]*LastOut[1][s2]
				+L_a[k-1]-log(1+exp(L_a[k-1]));
			for (s1=0;s1<nstate;s1++)
				sum+=exp(gamma[s1]+Alpha[k-1][s1]);
			if (sum<1E-300)
				Alpha[k][s2]=-INIFINITY_;
			else
				Alpha[k][s2]=log(sum);
		}
		// normalization,prevent overflow
		max_branch[k]=Alpha[k][0];
		for (s2=1;s2<nstate;s2++)
			if (Alpha[k][s2]>max_branch[k])
				max_branch[k]=Alpha[k][s2];
		for (s2=0;s2<nstate;s2++)
			Alpha[k][s2]=Alpha[k][s2]-max_branch[k];
	}

	// backward recursion,compute Beta
	for (k=L_TOTAL-1;k>0;k--)
	{
		for (s1=0;s1<nstate;s1++)
		{
			for (s2=0;s2<nstate;s2++)	// initialize metric
				gamma[s2]=-INIFINITY_;
			gamma[NextState[0][s1]]=-msg[k]+parity[k]*NextOut[0][s1]
				-log(1+exp(L_a[k]));	// bit0 
			gamma[NextState[1][s1]]=msg[k]+parity[k]*NextOut[1][s1]
				+L_a[k]-log(1+exp(L_a[k]));	// bit1 
			sum=0.0;
			for (s2=0;s2<nstate;s2++)
				sum+=exp(gamma[s2]+Beta[k+1][s2]);
			if (sum<1E-300)
				Beta[k][s1]=-INIFINITY_;
			else
				Beta[k][s1]=log(sum);
		}
		for (s1=0;s1<nstate;s1++)
			Beta[k][s1]=Beta[k][s1]-max_branch[k];
	}

	// forward again,computer LLRs
	for (k=0;k<L_TOTAL;k++)
	{
		sum0=0;sum1=0;
		for (s2=0;s2<nstate;s2++)
		{
			gamma[LastState[0][s2]]=-msg[k]+parity[k]*LastOut[0][s2]
				-log(1+exp(L_a[k]));
			gamma[LastState[1][s2]]=msg[k]+parity[k]*LastOut[1][s2]
				+L_a[k]-log(1+exp(L_a[k]));
			sum0+=exp(gamma[LastState[0][s2]]+Alpha[k][LastState[0][s2]]+Beta[k+1][s2]);
			sum1+=exp(gamma[LastState[1][s2]]+Alpha[k][LastState[1][s2]]+Beta[k+1][s2]);
		}
		L_all[k]=log(sum1)-log(sum0);
	}

	delete [] Beta;
	delete [] Alpha;
	delete [] max_branch;
}

//////////////////////////////////////////////////////////////////////
// universal turbo code decoder
// stream -- received stream, scaled according to estimated CSI
// msghat -- {0,1},deocder output, estimated information bits 
// puncture -- if punctured, NOT tested yet
// niter -- max iteration number
// stoprule -- stop rule, 0=max iter;1=CRC;2=SDR;
//////////////////////////////////////////////////////////////////////
UINT decode(double *stream, BYTE *msghat,BYTE *m,int *bits_err,int *frame_err,BOOL puncture, BOOL usemap,INT niter, INT stoprule)
{
	double msg[L_TOTAL],imsg[L_TOTAL];
	double parity[2][L_TOTAL];
	double L_all[L_TOTAL];
	double L_e[L_TOTAL],L_a[L_TOTAL];
	INT i,iter;
	bool f_err;

	// demultiplex
	for (i=0;i<L_TOTAL;i++)
	{
		if (puncture){	// punctured rate=1/2
			msg[i]=stream[2*i];
			parity[i%2][i]=stream[i*2+1];
		}else{			// unpunctured rate=1/3
			msg[i]=stream[3*i];
			parity[0][i]=stream[3*i+1];
			parity[1][i]=stream[3*i+2];
		}
	}

	// interleave system bit for decoder 2#
	for (i=0;i<L_TOTAL;i++)
		imsg[i]=msg[m_Inter_table[i]];

	// initialize extrinsic information
	for (i=0;i<L_TOTAL;i++)
		L_e[i]=0;
	
	// start decode iteration
	for (iter=0;iter<niter;iter++)
	{
		f_err=false;
		// decoder 1#
		for (i=0;i<L_TOTAL;i++)	//deinterleave a prior information
			L_a[m_Inter_table[i]]=L_e[i];
		if (usemap)
			logmap(msg,parity[0],L_a,L_all,true);
		// compute extrinsic information
		for (i=0;i<L_TOTAL;i++)
			L_e[i]=L_all[i]-2*msg[i]-L_a[i];

		// decoder 2#
		for (i=0;i<L_TOTAL;i++)	// interleave a prior information
			L_a[i]=L_e[m_Inter_table[i]];
		if (usemap)
			logmap(imsg,parity[1],L_a,L_all,false);

		// compute extrinsic information
		for (i=0;i<L_TOTAL;i++)
			L_e[i]=L_all[i]-2*imsg[i]-L_a[i];

		// estimate information bits
		for (i=0;i<L_TOTAL;i++)
			if(L_all[i]>0)
				msghat[m_Inter_table[i]]=1;
			else
				msghat[m_Inter_table[i]]=0;
		// count errors
		for (i=0;i<L_TOTAL-M;i++)
		{
			if (msghat[i]!=m[i])
			{
				bits_err[iter]++;
				f_err=true;
			}
		}
		if (f_err)
			frame_err[iter]++;
		// check stop rule
		// expected,:)
	}
	return iter;
}



int main(int argc, char* argv[])
{
    initRandom(seed);

	BYTE * m;
	BYTE * x;
	double * y;
	BYTE * mhat;

	int frame;
	int bits_all,bits_err[MAXITER],frame_err[MAXITER];
	double Ber,Fer;
	double Eb_No_dB,No;
	bool f_err;
	FILE * fp;
	int i;

	m = new BYTE[L_TOTAL];
	x = new BYTE[L_ALL];
	y = new double[L_ALL];
	mhat = new BYTE[L_TOTAL];


	if ((fp=fopen("LOGMAP.txt","w"))==NULL)
	{
		printf("\nFile open error!");
		return 1;
	}
	fprintf(fp,"Turbo Code Performance Simulation in Gaussian Channel\n");
	fprintf(fp,"rate=1/3,LogMap used\n");


	//initInter_table();
	init_Block_interleave_table();	// block interleave

	for (Eb_No_dB=-3.0;Eb_No_dB<5.0;Eb_No_dB+=0.5)
	{
		// dB = 10*log(Eb/No) where Eb is 1
		No = 1/pow(10.0,Eb_No_dB/10.0);
		
		bits_all=0;
		for (i=0;i<MAXITER;i++)
		{
			bits_err[i]=0;
			frame_err[i]=0;
		}

		for (frame=0;frame<FRAME_NUM;frame++,bits_all+=L_TOTAL-M)
		{
			f_err=false;
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
			// decode
			decode(y,mhat,m,bits_err,frame_err,false,true,MAXITER,0);
			// count errors

		}	// end for --- one frame

		printf("-------------------------\n");
		printf("Eb/No=%fdB:\n",Eb_No_dB);
		printf("-------------------------\n");
		fprintf(fp,"-------------------------\n");
		fprintf(fp,"Eb/No=%fdB:\n",Eb_No_dB);
		fprintf(fp,"-------------------------\n");

		for (i=0;i<MAXITER;i++)
		{
			Ber=(double)bits_err[i]/(double)bits_all;
			Fer=(double)frame_err[i]/(double)FRAME_NUM;
			printf("Iteration:%d\n",i+1);
			printf("---Ber=%f\n---Fer=%f\n",Ber,Fer);
			fprintf(fp,"Iteration:%d\n",i);
			fprintf(fp,"---Ber=%f\n---Fer=%f\n",Ber,Fer);
		}
	}	// end for ---one Eb/No,assigned total bits

	fclose(fp);

	delete m;
	delete x;
	delete y;
	delete mhat;

	return 0;
}

