#include "log_map.h"
#include <fstream>
#include <iostream>
using namespace std;

extern int length_after_code;
extern int source_length;
extern int MODULATION;
extern int SYMBOL_NUM;
extern int f1,f2;

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


