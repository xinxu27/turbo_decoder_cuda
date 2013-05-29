#include "main.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
using namespace std;

/*==================================================*/
/* ������� */
#define TERMINATED          1               /* 0����β��1��β */

#define TYPE_DECODER		1				/* ����������:1-LogMAP����
														  2-MAX-LogMAP����
														  3-SOVA����
                                                          4-const_LogMAP */
#define N_ITERATION			8				/* ����������� */
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
/* Log-MAP�㷨�õ��Ĳ��ұ� */
const double lookup_index_Log_MAP[16] = {0.0, 0.08824, 0.19587, 0.31026, 0.43275, 0.56508,
								0.70963, 0.86972, 1.0502, 1.2587, 1.5078, 1.8212,
								2.2522, 2.9706, 3.6764, 4.3758};
const double lookup_table_Log_MAP[16] = {0.69315, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35,
								0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.025, 0.0125};

/* Log-MAP�㷨�õ��Ĳ��ұ� */
const long lookup_index_Log_MAP_fix[9] = {0,1,2,3,4,5,6,7,8};
const long lookup_table_Log_MAP_fix[9] = {3,2,2,2,1,1,1,1,1};

//const long lookup_index_Log_MAP_fix[9] = {0,1,2,3,4,5,6,7,8};
//const long lookup_table_Log_MAP_fix[9] = {0,0,0,0,0,0,0,0,0};
/*==================================================*/
/* ����������صĲ��� */
int M_num_reg = COLUMN_OF_G-1;		/* �Ĵ����� */
int n_states = 8;						/* ״̬��:2��M_num_reg���� */
/*==================================================*/

//int M_num_reg;
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

double _qpsk_map_i[4]=
{
	0.7071,0.7071,-0.7071,-0.7071
};

double _qpsk_map_q[4]=
{
	0.7071,-0.7071,0.7071,-0.7071
};
double _8psk_map_i[8]=
{
	-0.7071,-1,0,0.7071,0,-0.7071,0.7071,1
};

double _8psk_map_q[8]=
{
	0.7071,0,1,0.7071,-1,-0.7071,-0.7071,0
};
double _16QAM_map_i[16]=
{
	-0.948683,-0.948683,-0.948683,-0.948683,
	-0.316228,-0.316228,-0.316228,-0.316228,
	0.948683,0.948683,0.948683,0.948683,
	0.316228,0.316228,0.316228,0.316228,
};

double _16QAM_map_q[16]=
{
	-0.948683,-0.316228,0.948683,0.316228,
	-0.948683,-0.316228,0.948683,0.316228,
	-0.948683,-0.316228,0.948683,0.316228,
	-0.948683,-0.316228,0.948683,0.316228,

};
double _64QAM_map_i[64]=
	{
	   0.4629,0.4629,0.4629,0.4629,0.4629,0.4629,0.4629,0.4629,0.1543,0.1543,0.1543,0.1543,0.1543,0.1543,0.1543,0.1543,
	   0.7615,0.7615,0.7615,0.7615,0.7615,0.7615,0.7615,0.7615,1.0801,1.0801,1.0801,1.0801,1.0801,1.0801,1.0801,1.0801,
	   -0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.4629,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,-0.1543,
	   -0.7615,-0.7615,-0.7615,-0.7615,-0.7615,-0.7615,-0.7615,-0.7615,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801,-1.0801
	
	
	};


double _64QAM_map_q[64]=
	{
	    0.4629,0.1543,0.7615,1.0801,-0.4629,-0.1543,-0.7615,-1.0801,0.4629,0.1543,0.7615,1.0801,-0.4629,-0.1543,-0.7615,-1.0801,
	    0.4629,0.1543,0.7615,1.0801,-0.4629,-0.1543,-0.7615,-1.0801,0.4629,0.1543,0.7615,1.0801,-0.4629,-0.1543,-0.7615,-1.0801,
	    0.4629,0.1543,0.7615,1.0801,-0.4629,-0.1543,-0.7615,-1.0801,0.4629,0.1543,0.7615,1.0801,-0.4629,-0.1543,-0.7615,-1.0801,
	    0.4629,0.1543,0.7615,1.0801,-0.4629,-0.1543,-0.7615,-1.0801,0.4629,0.1543,0.7615,1.0801,-0.4629,-0.1543,-0.7615,-1.0801
	
	
	};
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




void _qpsk_module(int * a,double * outi,double * outq,int N)
{
	int n=N/2;

	int temp,i;

	for (i=0;i<n;i++)
	{    
		temp=(*(a+i*2))*2+(*(a+i*2+1));
			
		*(outi+i)=*(_qpsk_map_i+temp);
			
		*(outq+i)=*(_qpsk_map_q+temp);
	}
 
}




void _8psk_module(int * a,double * outi,double * outq,int N)
{
	int n=N/3;

	int temp,i;
	
	for (i=0;i<n;i++)
	{    
		temp=((*(a+i*3+2))*4+(*(a+i*3+1))*2+*(a+i*3));
			
		*(outi+i)=*(_8psk_map_i+temp);
			
		*(outq+i)=*(_8psk_map_q+temp);
	} 
}



void _16qam_module(int * a,double * outi,double * outq,int N)
{
	int n=N/4;
	int temp,i;
	for (i=0;i<n;i++)
	{    
		temp=(*(a+i*4))*8+(*(a+i*4+1))*4+(*(a+i*4+2))*2+*(a+i*4+3);
			
		*(outi+i)=*(_16QAM_map_i+temp);
			
		*(outq+i)=*(_16QAM_map_q+temp);
	}
}
void _64qam_module(int * a,double * outi,double * outq,int N)
{
	int n=N/6;
	int temp,i;
	for (i=0;i<n;i++)
	{    
	     
		temp=(*(a+i*6))*32+(*(a+i*6+1))*16+(*(a+i*6+2))*8+*(a+i*6+3)*4+*(a+i*6+4)*2+*(a+i*6+5);
			
		*(outi+i)=*(_64QAM_map_i+temp);
			
		*(outq+i)=*(_64QAM_map_q+temp);
	}
}



void module(int * a,double * outi,double * outq,int N, int modu_index)
{
	switch (modu_index)
	{
	    case 1: _bpsk_module(a,outi,outq,N);break;
		case 2: _qpsk_module(a,outi,outq,N);break;
	    case 3: _8psk_module(a,outi,outq,N);break;
		case 4: _16qam_module(a,outi,outq,N);break;
		case 6: _64qam_module(a,outi,outq,N);break;
	}

}


void _bpsk_demodule(double *symbol_i, double *symbol_q, int symbol_len,
					double *out,double Kf)
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



void _qpsk_demodule(double *symbol_i, double *symbol_q, int symbol_len,
					double *out,double Kf)
{
	int i,j;
	double min_sqr_dis1,min_sqr_dis2,sqr_dis1,sqr_dis2;

	for(i=0;i<symbol_len;i++)
	{
	/*the minimum square distance has nothing with the orders of calculation of sqr_dis1 and sqr_dis2*/

		min_sqr_dis1=0x7fffffffffff,min_sqr_dis2=0x7fffffffffff;

        for(j=0;j<4;j++)
		{
			if(j&2)
			{
			/*calculate Min{dj} where j in Sj*/
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_qpsk_map_i[j],(double)_qpsk_map_q[j]);
	
			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_qpsk_map_i[j],(double)_qpsk_map_q[j]);

			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}
		}

		out[2*i]=-Kf*(min_sqr_dis1-min_sqr_dis2);

		min_sqr_dis1=0x7fffffffffff,min_sqr_dis2=0x7fffffff;

		for(j=0;j<4;j++)
		{
			if(j&1)
			{
			/*calculate Min{dj} where j in Sj*/
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_qpsk_map_i[j],(double)_qpsk_map_q[j]);
	
			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_qpsk_map_i[j],(double)_qpsk_map_q[j]);

			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}

		out[2*i+1]=-Kf*(min_sqr_dis1-min_sqr_dis2);
	}
}

void _8psk_demodule(double *symbol_i, double *symbol_q, int symbol_len,
					double *out,double Kf)
{
	int i,j;
	double min_sqr_dis1,min_sqr_dis2,sqr_dis1,sqr_dis2;

	for(i=0;i<symbol_len;i++)
	{
	/*the minmun square distance has nothing with the orders of calculation of sqr_dis1 and sqr_dis2*/
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;

		for(j=0;j<8;j++)
		{
			if(j&4)
			{
			/*calculate Min{dj} where j in Sj*/
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_8psk_map_i[j],(double)_8psk_map_q[j]);
	
			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_8psk_map_i[j],(double)_8psk_map_q[j]);

			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}

		out[3*i+2]=-Kf*(min_sqr_dis1-min_sqr_dis2);
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;

		for(j=0;j<8;j++)
		{
			if(j&2)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_8psk_map_i[j],(double)_8psk_map_q[j]);	

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_8psk_map_i[j],(double)_8psk_map_q[j]);	

			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[3*i+1]=-Kf*(min_sqr_dis1-min_sqr_dis2);
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;
		
		for(j=0;j<8;j++)
		{
			if(j&1)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_8psk_map_i[j],(double)_8psk_map_q[j]);					

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_8psk_map_i[j],(double)_8psk_map_q[j]);	
			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[3*i ]=-Kf*(min_sqr_dis1-min_sqr_dis2);
	}
}

void _16qam_demodule(double *symbol_i, double *symbol_q, int symbol_len,
					double* out,double Kf)
{
    int i,j;
	double min_sqr_dis1,min_sqr_dis2,sqr_dis1,sqr_dis2;

	for(i=0;i<symbol_len;i++)
	{
	/*the minmun square distance has nothing with the orders of calculation of sqr_dis1 and sqr_dis2*/
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;

		for(j=0;j<16;j++)
		{
			if(j&8)
			{
			/*calculate Min{dj} where j in Sj*/
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_16QAM_map_i[j],(double)_16QAM_map_q[j]);
	
			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_16QAM_map_i[j],(double)_16QAM_map_q[j]);

			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}

		out[4*i]=-Kf*(min_sqr_dis1-min_sqr_dis2);
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;

		for(j=0;j<16;j++)
		{
			if(j&4)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_16QAM_map_i[j],(double)_16QAM_map_q[j]);	

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_16QAM_map_i[j],(double)_16QAM_map_q[j]);	

			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[4*i+1]=-Kf*(min_sqr_dis1-min_sqr_dis2);
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;
		
		for(j=0;j<16;j++)
		{
			if(j&2)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_16QAM_map_i[j],(double)_16QAM_map_q[j]);					

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_16QAM_map_i[j],(double)_16QAM_map_q[j]);	
			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[4*i+2]=-Kf*(min_sqr_dis1-min_sqr_dis2);
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;

		for(j=0;j<16;j++)
		{
			if(j&1)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_16QAM_map_i[j],(double)_16QAM_map_q[j]);					

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_16QAM_map_i[j],(double)_16QAM_map_q[j]);	
			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[4*i+3]=-Kf*(min_sqr_dis1-min_sqr_dis2);
	}
}   


void _64qam_demodule(double *symbol_i, double *symbol_q, int symbol_len,
					double* out,double Kf)
{
	int i,j;
	double min_sqr_dis1,min_sqr_dis2,sqr_dis1,sqr_dis2;

	for(i=0;i<symbol_len;i++)
	{
	/*the minmun square distance has nothing with the orders of calculation of sqr_dis1 and sqr_dis2*/
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;

		for(j=0;j<64;j++)
		{
			if(j&32)
			{
			/*calculate Min{dj} where j in Sj*/
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);
	
			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);

			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}

		out[6*i]=-Kf*(min_sqr_dis1-min_sqr_dis2);
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;

		for(j=0;j<64;j++)
		{
			if(j&16)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);	

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);	

			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[6*i+1]=-Kf*(min_sqr_dis1-min_sqr_dis2);
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;
		
		for(j=0;j<64;j++)
		{
			if(j&8)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);					

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);	
			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[6*i+2]=-Kf*(min_sqr_dis1-min_sqr_dis2);
		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;

		for(j=0;j<64;j++)
		{
			if(j&4)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);					

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);	
			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[6*i+3]=-Kf*(min_sqr_dis1-min_sqr_dis2);

		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;
		
		for(j=0;j<64;j++)
		{
			if(j&2)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);					

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);	
			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[6*i+4]=-Kf*(min_sqr_dis1-min_sqr_dis2);

		min_sqr_dis1=0x7fffffff,min_sqr_dis2=0x7fffffff;
		
		for(j=0;j<64;j++)
		{
			if(j&1)
			{
				/*calculate Min{dj} where j in Sj*/		
			sqr_dis1=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);					

			if(sqr_dis1<min_sqr_dis1)
				min_sqr_dis1=sqr_dis1;
			}

			else
			{
				/*calculate Min{dj} where j in complement of Sj*/	
			sqr_dis2=calculate_sqr_dis(symbol_i[i],symbol_q[i],(double)_64QAM_map_i[j],(double)_64QAM_map_q[j]);	
			if(sqr_dis2<min_sqr_dis2)
				min_sqr_dis2=sqr_dis2;

			}


		}
		out[6*i+5]=-Kf*(min_sqr_dis1-min_sqr_dis2);
	
	


	}


}


void demodule(double *symbol_i, double *symbol_q, int symbol_len,double* out,double Kf,int modu_index)
{
    switch(modu_index)
	{
	    case 1: _bpsk_demodule(symbol_i, symbol_q,symbol_len,out,Kf);break;
		case 2: _qpsk_demodule(symbol_i, symbol_q,symbol_len,out,Kf);break;
	    case 3: _8psk_demodule(symbol_i, symbol_q,symbol_len,out,Kf);break;
		case 4: _16qam_demodule(symbol_i, symbol_q,symbol_len,out,Kf);break;
		case 6: _64qam_demodule(symbol_i, symbol_q,symbol_len,out,Kf);break;
		default:printf("please input the number between 1 and 6.");
	}
}





int main()
{
/*void main(int argc, char *argv[])
{
          int RB=atoi(argv[1]);
	int MCS=atoi(argv[2]);
	MODULATION = atoi(argv[3]);	//6;			//���ƽ�����1,2,3,4,6
	length_after_code = atoi(argv[4]);//SYMBOL_NUM*MODULATION;	//������׺�ĳ���	
	SYMBOL_NUM=length_after_code/MODULATION;//  atoi(argv[4]);//672

	source_length =atoi(argv[5]);//2688;			//������Դ����
	f1 =atoi(argv[6]);// 127;
	f2 =atoi(argv[7]);//504; 					//���볤�йصĽ�֯������

	double EbN0start = atof(argv[8]);//0;		//������ʼ�����
	double EbN0end = atof(argv[9]);//4;		//��������ֹ�����
	double EbN0step =atof(argv[10]);;		//��������Ȳ���
 */
    
	
	//int RB=1;
	//int MCS=atoi(argv[2]);
	MODULATION = 1;			//���ƽ�����1,2,3,4,6
	source_length = 6144;			//������Դ����
	length_after_code = 6144*3+12;//SYMBOL_NUM*MODULATION;	//������׺�ĳ���	
	//SYMBOL_NUM=length_after_code/MODULATION;//  atoi(argv[4]);//672
	SYMBOL_NUM=length_after_code;

	
	f1 =263;
	f2 =480; 					//���볤�йصĽ�֯������

	double EbN0start = 0;//0;		//������ʼ�����
	double EbN0end = 2;//4;		//��������ֹ�����
	double EbN0step =0.2;		//��������Ȳ���
	int FrameNum=100;	

	int floatorfix = 0;//1-fix,0-float
		
	//rate=(double)source_length/(double)length_after_code;
	double rate = (double)source_length/(double)(SYMBOL_NUM);


	int *source = NULL;

	int *coded_source = NULL;
	int *punced_source = NULL;
// int out_len;
	double *modulated_source_i, *modulated_source_q;

	double *after_channel_i,*after_channel_q;
	double *received_punced_source = NULL;
	double *flow_for_decode = NULL;
	int *flow_decoded = NULL;

	double EbN0dB,sigma;
	int nf, i1,snr_num,result_num=0;
	double *err_bit_rate,*err_block_rate;

	//int N_ITERATION=8;
	int temp[8],err_bit_num[8], err_block_num[8];
	int i2;

	FILE *fp;   

	
	if ((fp=fopen("result.txt","a+"))==NULL)
	{
		printf("cannot open file\n");
	    return 0;
	}

	//fprintf(fp,"\n%s","MCS = ");
	//fprintf(fp,"%d", MCS);

	fprintf(fp,"\n%s","MODULATION = ");
	fprintf(fp,"%d", MODULATION);

	//fprintf(fp,"\n%s","RB = ");
	//fprintf(fp,"%d", RB);

	fprintf(fp,"\n%s","source_length = ");
	fprintf(fp,"%d", source_length);

	//fprintf(fp,"\n%s","transport_bit_length = ");
	//fprintf(fp,"%d", length_after_code);
	
	fprintf(fp,"\n%s","rate_coding = ");
	fprintf(fp,"%f",  rate);
	

	fprintf(fp,"\n%s", "SNR_begin=");
	fprintf(fp, "%f", EbN0start);

	fprintf(fp,"\n%s", "SNR_step=");
	fprintf(fp,"%f", EbN0step);



	TurboCodingInit();
			
	if ((source=(int *)malloc(source_length*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of source \n");
	  exit(1);  
	}
	if ((coded_source=(int *)malloc((3*source_length+4*M_num_reg)*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of coded_source \n");
	  exit(1);  
	}
	/*if ((punced_source=(int *)malloc((int)(SYMBOL_NUM*MODULATION)*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of punced_source \n");
	  exit(1);  
	}*/
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
	/*if ((received_punced_source=(double *)malloc((int)(SYMBOL_NUM*MODULATION)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of received_punced_source \n");
	  exit(1);  
	}*/
	if ((flow_for_decode=(double *)malloc((3*source_length+4*M_num_reg)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of flow_for_decode \n");
	  exit(1);  
	}
	if ((flow_decoded=(int *)malloc(N_ITERATION*source_length*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of flow_decoded\n");
	  exit(1);  
	}
	snr_num=(int)((EbN0end-EbN0start)/EbN0step+1);
	if ((err_bit_rate=(double*)malloc(N_ITERATION*snr_num*sizeof(double)))==NULL)
	{
		printf("\n fail to allocate memory of err_bit_rate \n");
		exit(1);
	}
	if ((err_block_rate=(double*)malloc(N_ITERATION*snr_num*sizeof(double)))==NULL)
	{
		printf("\n fail to allocate memory of err_block_rate \n");
		exit(1);
	}
	
	srand((unsigned)time(NULL));

	for (EbN0dB=EbN0start; EbN0dB<=EbN0end; EbN0dB+=EbN0step)
	{
		sigma = pow(10,-EbN0dB/20)*sqrt(0.5/(rate*MODULATION));
		//sigma = pow(10,-EbN0dB/20)*sqrt(0.5);
		for (i2=0;i2<N_ITERATION;i2++)
		{
			temp[i2]=0;err_bit_num[i2]=0; err_block_num[i2]=0;
		}

		for (nf=0; nf<FrameNum; nf++)
		{
			for(i1=0; i1<source_length; i1+=1)
			{
				*(source+i1)=rand()%2;
				//*(source+i1)=(int)(rand()/16384);
			}




            TurboEnCoding(source, coded_source, source_length);

/*******************************************************************************/
			//rate match
			//rate_match(coded_source,3*source_length+4*M_num_reg,punced_source,length_after_code);
			module(coded_source,modulated_source_i,modulated_source_q,SYMBOL_NUM*MODULATION,MODULATION);

			AWGN(modulated_source_i, after_channel_i, sigma, SYMBOL_NUM);
			AWGN(modulated_source_q, after_channel_q, sigma, SYMBOL_NUM);

			demodule(after_channel_i, after_channel_q, SYMBOL_NUM,flow_for_decode,1/(2*pow(sigma,2)),MODULATION);
			//de rate match
			//de_rate_match(received_punced_source,flow_for_decode,length_after_code,3*source_length+4*M_num_reg);
/*******************************************************************************/

		/*ofstream flow_for_decode_f("flow_for_decode.txt");
		for(int j=0;j<3*source_length+4*M_num_reg;j++)
		{
			flow_for_decode_f << flow_for_decode[j] << endl;
		}*/

		ifstream flow_for_decode_reader("flow_for_decode.txt");
		for (int j=0; j<3*source_length+4*M_num_reg; j++) 
		{
			flow_for_decode_reader >> flow_for_decode[j];
		}



		TurboDecoding(flow_for_decode, flow_decoded, 3*source_length+4*M_num_reg);
			

		for(i2=0;i2<N_ITERATION;i2++)
		{
				temp[i2] = err_bit_num[i2];
			
				for (i1=0; i1<source_length; i1++)
				{
					if (*(source+i1) != *(flow_decoded+i2*source_length+i1))
					{
						err_bit_num[i2] = err_bit_num[i2]+1;
					}
				}
				if(temp[i2]!=err_bit_num[i2])
					err_block_num[i2]++;
			}
			
				if(err_block_num[N_ITERATION-1]>=50)//����1000���飬����
				{
					nf++;
					break;
				}
		}
		for(i2=0;i2<N_ITERATION;i2++)		
		{
			err_bit_rate[result_num*N_ITERATION+i2] = (double)err_bit_num[i2]/(nf*source_length);
			err_block_rate[result_num*N_ITERATION+i2] = (double)err_block_num[i2]/(double)(nf);
		
/*-----------------------------------------------------------------*/
			printf("block error rate: ");
			printf("%.10f ",err_block_rate[result_num*N_ITERATION+i2]);
			printf("bit error rate: ");
			printf("%.10f ",err_bit_rate[result_num*N_ITERATION+i2]);
			printf("\n");
		
		}
	//	if (err_block_rate[result_num]<=0.01)
	//	{
	//		result_num++;
	//		break;//�����С��0.01������
	//	}
		printf("\n");
		result_num++;

	}
    fprintf(fp,"%s ","N_ITERATION=1");
	if (EbN0dB>EbN0end)
	{
		fprintf(fp,"\n%s", "SNR_end=");
		fprintf(fp,"%f", EbN0end);
	}
	else
	{
		fprintf(fp,"\n%s", "SNR_end=");
		fprintf(fp,"%f", EbN0dB);
	}
	double factor;
	factor=10*log(rate*MODULATION)/log(10);
	fprintf(fp,"\n%s\n", "Es/N0:");
	for (i1=0;i1<result_num;i1++)
	{
		fprintf(fp,"%f ",(EbN0start+EbN0step*i1)+factor);
	}

	fprintf(fp,"\n%s\n", "Ber:");
	
	for(i2=0;i2<N_ITERATION;i2++)
	{
		for(i1=0;i1<result_num;i1++)
		{
			fprintf(fp," %.10f ",err_bit_rate[i1*N_ITERATION+i2]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n%s\n", "Bler:");
	for(i2=0;i2<N_ITERATION;i2++)
	{
		for(i1=0;i1<result_num;i1++)
		{
			fprintf(fp,"%.10f ",err_block_rate[i1*N_ITERATION+i2]);
		}
		fprintf(fp,"\n");
		
	}
	fprintf(fp,"\n%s\n", "throughput:");
	for (i1=0;i1<result_num;i1++)
	{
		fprintf(fp,"%.10f ",(1-err_block_rate[i1])*rate*MODULATION);
	}


fprintf(fp,"----------------------------------------------------------");

	fclose(fp);
/*-----------------------------------------------------------------*/	
/*-----------------------------------------------------------------*/	
	TurboCodingRelease();
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/	

	free(source);
	free(coded_source);
	//free(punced_source);
	free(modulated_source_i);
	free(modulated_source_q);
	free(after_channel_i);
	free(after_channel_q);
	//free(received_punced_source);
	free(flow_for_decode);
	free(flow_decoded);
}


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
double E_algorithm(double x, double y)
{
	double temp = (y-x)>0? (y-x):(x-y);
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
		temp = (double)lookup_table_Log_MAP[i-1];
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
double E_algorithm_seq(double *data_seq, int length)
{
	int i;			/* ѭ������ */
	double temp;

	/* ÿ��������E�㷨,������һ������E�㷨 */
	temp = E_algorithm(*(data_seq+0), *(data_seq+1));
	for (i=2; i<length; i++)
	{
		temp = E_algorithm(temp, *(data_seq+i));
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

/*---------------------------------------------------------------
����:
	void Log_MAP_decoder(double *recs_turbo, double *La_turbo, int terminated, 
						double *LLR_all_turbo, int len_total)
����:
	Log_MAP������.
����:
	�������:
		rec_turbo -��Lc��Ľ���������ַ.
		La_turbo - ����Ϣ��
		terminated - �Ƿ��β.
		len_total - �������г���.
	�������:
		LLR_all_turbo - ��Ȼ��.
����ֵ:
	�ޣ�
---------------------------------------------------------------*/
void Log_MAP_decoder(double *recs_turbo, double *La_turbo, int terminated, double *LLR_all_turbo, int len_total)
{
	int i, j;			/* ѭ������ */

	double *alpha_Log, *beta_Log;		/* n_states*(len_total+1) */
	double *gama_Log;					/* n_states*len_total*2 */

	double *tempmax;					/* len_total+1 */
	double *temp0, *temp1;				/* n_states */
	double tempx, tempy;

	if ((alpha_Log=(double *)malloc(n_states*(len_total+1)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of alpha_Log \n");
	  exit(1);  
	}
	if ((beta_Log=(double *)malloc(n_states*(len_total+1)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of beta_Log \n");
	  exit(1);  
	}
	if ((gama_Log=(double *)malloc(n_states*len_total*2*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of gama_Log \n");
	  exit(1);  
	}

	if ((tempmax=(double *)malloc((len_total+1)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of tempmax \n");
	  exit(1);  
	}

	if ((temp0=(double *)malloc(n_states*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of temp0 \n");
	  exit(1);  
	}
	if ((temp1=(double *)malloc(n_states*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of temp1 \n");
	  exit(1);  
	}
	
	/*===== ��ʼ��alpha��beta =====*/
	*(alpha_Log+0) = 0;
	*(beta_Log+len_total) = 0;

	for (i=1; i<n_states; i++)
	{
		*(alpha_Log+i*(len_total+1)+0) = (double)-INFTY;

		/* ��β */
		if (terminated)
		{
			*(beta_Log+i*(len_total+1)+len_total) = (double)-INFTY;
		}
		/* ����β */
		else
		{
			*(beta_Log+i*(len_total+1)+len_total) = 0;
		}
	}
	/*======== ����Gama ========*/
	for (i=0; i<len_total; i++)		/* �����ѭ�� */
	{
		for (j=0; j<n_states; j++)	/* ��״̬ѭ�� */
		{
			/* ������gama����ʽ */
			*(gama_Log+j*len_total*2+i*2+0) 
				= -*(recs_turbo+2*i) + *(recs_turbo+2*i+1)*(*(turbo_trellis.mx_nextout+j*4+1)) - *(La_turbo+i)/2;
			*(gama_Log+j*len_total*2+i*2+1)
				= *(recs_turbo+2*i) + *(recs_turbo+2*i+1)*(*(turbo_trellis.mx_nextout+j*4+3)) + *(La_turbo+i)/2;
		}	/* ��״̬ѭ������ */
	}		/* �����ѭ������ */

	/*======== ǰ�����alpha ========*/
	for (i=1; i<len_total+1; i++)	/* �����ѭ�� */
	{
		for (j=0; j<n_states; j++)	/* ��״̬ѭ�� */
		{
			tempx = *(gama_Log+(*(turbo_trellis.mx_laststat+j*2+0))*len_total*2+(i-1)*2+0)
					+ *(alpha_Log+(*(turbo_trellis.mx_laststat+j*2+0))*(len_total+1)+i-1);
			tempy = *(gama_Log+(*(turbo_trellis.mx_laststat+j*2+1))*len_total*2+(i-1)*2+1)
					+ *(alpha_Log+(*(turbo_trellis.mx_laststat+j*2+1))*(len_total+1)+i-1);
			*(alpha_Log+j*(len_total+1)+i) = E_algorithm(tempx, tempy);
		}	/* ��״̬ѭ������ */

		/* ����tempmax[i],���ڹ�һ��alpha��beta */
		for (j=0; j<n_states; j++)
		{
			if (*(tempmax+i) < *(alpha_Log+j*(len_total+1)+i))
			{
				*(tempmax+i) = *(alpha_Log+j*(len_total+1)+i);
			}
		}

		/* ��һ��alpha */
		for(j=0; j<n_states; j++)
		{
			*(alpha_Log+j*(len_total+1)+i) = *(alpha_Log+j*(len_total+1)+i) - *(tempmax+i);
		}

	}	/* �����ѭ������ */

	/*======== �������beta ========*/
	for (i=len_total-1; i>=0; i--)	/* �����ѭ�� */
	{
		for (j=0; j<n_states; j++)	/* ��״̬ѭ�� */
		{
			tempx = *(gama_Log+j*len_total*2+i*2+0) 
					+ *(beta_Log+(*(turbo_trellis.mx_nextstat+j*2+0))*(len_total+1)+i+1);
			tempy = *(gama_Log+j*len_total*2+i*2+1) 
					+ *(beta_Log+(*(turbo_trellis.mx_nextstat+j*2+1))*(len_total+1)+i+1);

			*(beta_Log+j*(len_total+1)+i) = E_algorithm(tempx, tempy);
		}	/* ��״̬ѭ������ */

		/* ��һ��beta */
		for (j=0; j<n_states; j++)
		{
			*(beta_Log+j*(len_total+1)+i) = *(beta_Log+j*(len_total+1)+i) - *(tempmax+i+1);
		}
	}	/* �����ѭ������ */

	/*=== ������Ȼ��LLR ===*/
	for (i=0; i<len_total; i++)		/* �����ѭ�� */
	{
		for (j=0; j<n_states; j++)	/* ��״̬ѭ�� */
		{
			*(temp0+j) = *(gama_Log+(*(turbo_trellis.mx_laststat+j*2+0))*len_total*2+i*2+0) 
				+ *(alpha_Log+*(turbo_trellis.mx_laststat+j*2+0)*(len_total+1)+i)
				+ *(beta_Log+ j*(len_total+1)+i+1);

			*(temp1+j) = *(gama_Log+(*(turbo_trellis.mx_laststat+j*2+1))*len_total*2+i*2+1) 
				+ *(alpha_Log+*(turbo_trellis.mx_laststat+j*2+1)*(len_total+1)+i)
				+ *(beta_Log+j*(len_total+1)+i+1);
		}	/* ��״̬ѭ������ */

		/* ������Ȼ�� */
		*(LLR_all_turbo+i) = E_algorithm_seq(temp1, n_states) - E_algorithm_seq(temp0, n_states);
	}	/* �����ѭ������ */

	free(alpha_Log);
	free(beta_Log);
	free(gama_Log);
	free(tempmax);
	free(temp0);
	free(temp1);
}

/*---------------------------------------------------------------
����:
	void Log_MAP_decoder_fix(int *recs_turbo, long *La_turbo, int terminated, 
						long *LLR_all_turbo, int len_total)
����:
	Log_MAP������.
����:
	�������:
		rec_turbo -��Lc��Ľ���������ַ.
		La_turbo - ����Ϣ��
		terminated - �Ƿ��β.
		len_total - �������г���.
	�������:
		LLR_all_turbo - ��Ȼ��.
����ֵ:
	�ޣ�
---------------------------------------------------------------*/



/*---------------------------------------------------------------
����:
	void demultiplex(double *rec_turbo, int len_info, double *yk_turbo)
����:
	�⸴��.
����:
	�������:
		rec_turbo -��������������ַ.
		len_info - �����������г��ȣ�
	�������:
		yk_turbo - �⸴�ú�����������ַ.
����ֵ:
	�ޣ�
---------------------------------------------------------------*/
void demultiplex(double *rec_turbo, int len_info, double *yk_turbo)
{
	int i;			/* ѭ������ */

	int len_total = len_info+M_num_reg;		/* �ܳ��� */

	double *info2, *inted_info2;			/* ��Ϣλ�ͽ�֯�����Ϣλ */

	if ((info2=(double *)malloc(len_info*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of info2 \n");
	  exit(1);  
	}
	if ((inted_info2=(double *)malloc(len_info*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of inted_info2 \n");
	  exit(1);  
	}
	
	/* ����Ϣ���� */
	for(i=0; i<len_info; i++)
	{
		*(info2+i) = *(yk_turbo+2*i) = *(rec_turbo+3*i);
		*(yk_turbo+2*i+1) = *(rec_turbo+3*i+1);
		*(yk_turbo+2*len_total+2*i+1) = *(rec_turbo+3*i+2);
	}

	/* ��֯��Ϣ���� */
	randominterleaver_double(info2, inted_info2, index_randomintlvr, len_info);

	for (i=0; i<len_info; i++)
	{
		*(yk_turbo+2*len_total+2*i) = *(inted_info2+i);
	}

	/* ��β���� */
	for (i=0; i<2*M_num_reg; i++)
	{
		*(yk_turbo+2*len_info+i) = *(rec_turbo+3*len_info+i);
		*(yk_turbo+2*len_total+2*len_info+i) = *(rec_turbo+3*len_info+2*M_num_reg+i);
	}

	free(info2);
	free(inted_info2);
}




/*---------------------------------------------------------------
����:
	void TurboDecoding(double *flow_for_decode, int *flow_decoded,
						  int *flow_length, double EbN0dB)
����:
	Turbo��ҵ����Ϣ���뺯��.
����:
	�������: flow_for_decode - ����������ַ.
			  flow_length - �������г���.
			  EbN0dB - �����.
	�������: flow_decoded - ������������ַ.
����ֵ:
	��
---------------------------------------------------------------*/
void TurboDecoding(double *flow_for_decode, int *flow_decoded,int flow_length)
{
	int i;							/* ѭ������ */
	int length_info, length_total;	/* ��Ϣλ�����ܳ��� */
	int iteration;					/* ����ѭ������ */

	float *receive_punc = NULL;					/* �������� */ 
	double *yk_turbo = NULL;					/* ����ys��yp */
	double *La_turbo, *Le_turbo, *LLR_all_turbo;		/* ����Ϣ����Ȼ�� */

	int *tempout;


	/* ������Ϣλ���� */
	length_info = (flow_length-4*M_num_reg)/3;

	/* �ܳ��� */
	length_total = length_info+M_num_reg;

	/* �����ڴ� */
	if ((receive_punc=(float *)malloc((3*length_info+4*M_num_reg)*sizeof(float)))==NULL)
	{
		printf("\n fail to allocate memory of receive_punc \n");
		exit(1);  
	}

	if ((yk_turbo=(double *)malloc(4*length_total*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of yk_turbo \n");
	  exit(1);  
	}

	if ((La_turbo=(double *)malloc(length_total*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of La_turbo \n");
	  exit(1);  
	}
	if ((Le_turbo=(double *)malloc(length_total*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of Le_turbo \n");
	  exit(1);  
	}
	if ((LLR_all_turbo=(double *)malloc(length_total*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of LLR_all_turbo \n");
	  exit(1);  
	}

	if ((tempout=(int *)malloc(length_total*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of tempout \n");
	  exit(1);  
	}

/************************************************************************/
//csqu-070717
	for (i=0; i<flow_length; i++)
	{
		*(flow_for_decode+i) *= 0.5;
	}
/************************************************************************/	

	/* �⸴�� */
	demultiplex(flow_for_decode, length_info, yk_turbo);	/* �⸴�� */
	
	/* ��ʼ������Ϣ����Ȼ��Ϣ */
	for (i=0; i<length_total; i++)
	{
		*(La_turbo+i) = *(Le_turbo+i) = *(LLR_all_turbo+i) = 0;
	}

	for (iteration=0; iteration<N_ITERATION; iteration++)		/* ��ʼ���� */
	{
		/* ������1: */
		/* �⽻֯����������2������Ϣ */
		random_deinterlvr_double(La_turbo, Le_turbo, index_randomintlvr, length_info);
		
		/* ��β����0 */
		for (i=length_info; i<length_total; i++)
		{
			*(La_turbo+i) = 0;
		}

			/* Log-MAP�㷨 */
			Log_MAP_decoder(yk_turbo, La_turbo, TERMINATED, LLR_all_turbo, length_total);
		

		/* ��������Ϣ */
		for (i=0; i<length_total; i++)
		{
			//csqu-070717
			*(Le_turbo+i) = *(LLR_all_turbo+i) - *(La_turbo+i) - 2*(*(yk_turbo+2*i));
		}

		/* ������2: */
		/* ��֯����������1������Ϣ */
		randominterleaver_double(Le_turbo, La_turbo, index_randomintlvr, length_info);

		/* ��β����0 */
		for (i=length_info; i<length_total; i++)
		{
			*(La_turbo+i) = 0;
		}
						
			/* Log-MAP�㷨 */
			Log_MAP_decoder(yk_turbo+2*length_total, La_turbo, TERMINATED, LLR_all_turbo, length_total);
		

		/* ��������Ϣ */
		for(i=0; i<length_total; i++)
		{
			//csqu-070717
			*(Le_turbo+i) = *(LLR_all_turbo+i) - *(La_turbo+i) - 2*(*(yk_turbo+2*length_total+2*i));
		}
		/* ������2���� */
		decision(LLR_all_turbo, length_total, tempout);

	/* ���о����ؽ⽻֯ */
		random_deinterlvr_int(flow_decoded+length_info*iteration, tempout, index_randomintlvr, length_info);
	}
	
	/* ��ȫ����Ȼ���о� */

	


	free(receive_punc);
	free(yk_turbo);

	free(La_turbo);
	free(Le_turbo);
	free(LLR_all_turbo);

	free(tempout);
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
		}

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
        }/*�������ļ��޶����������Ӹ�˹�ֲ��������*/
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
