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
/* Log-MAP算法用到的查找表 */
const double lookup_index_Log_MAP[16] = {0.0, 0.08824, 0.19587, 0.31026, 0.43275, 0.56508,
								0.70963, 0.86972, 1.0502, 1.2587, 1.5078, 1.8212,
								2.2522, 2.9706, 3.6764, 4.3758};
const double lookup_table_Log_MAP[16] = {0.69315, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35,
								0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.025, 0.0125};

/* Log-MAP算法用到的查找表 */
const long lookup_index_Log_MAP_fix[9] = {0,1,2,3,4,5,6,7,8};
const long lookup_table_Log_MAP_fix[9] = {3,2,2,2,1,1,1,1,1};

//const long lookup_index_Log_MAP_fix[9] = {0,1,2,3,4,5,6,7,8};
//const long lookup_table_Log_MAP_fix[9] = {0,0,0,0,0,0,0,0,0};
/*==================================================*/
/* 与生成阵相关的参数 */
int M_num_reg = COLUMN_OF_G-1;		/* 寄存器数 */
int n_states = 8;						/* 状态数:2的M_num_reg次幂 */
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
函数:
	int gen_g_matrix(int k_column, int g_row1, int g_row2, int *mx_g_turbo)
介绍:
	得到生成阵.
参数:
	输入参数:
		k_column - 生成阵列数.
		g_row1 - 生成阵第一行.
		g_row2 - 生成阵第二行.
	输出参数:
		mx_g_turbo - 生成阵首址.
返回值:
	1 - 成功地得到生成阵.
	0 - 生成生成阵失败.
---------------------------------------------------------------*/
int gen_g_matrix(int k_column, int g_row1, int g_row2, int *mx_g_turbo)
{
	int i, position;		/* 循环变量 */
	int high_num, low_num;
	
	/* 第一行 */
	high_num = g_row1;		
	position = 1;			/* 第几个8进制数 */
	while (high_num>0)
	{
		low_num = high_num%10;	/* 得到这个8进制数 */
		if (low_num>7)			/* 判断是否为8进制数 */
		{
			return 0;
		}
		high_num = high_num/10;		/* 记下其余部分 */

		/* 将8进制数转为二进制并保存 */
		for (i=k_column-(position-1)*3-1; i>=0 && i>=k_column-position*3; i--)
		{
			*(mx_g_turbo+i) = low_num%2;
			low_num = low_num/2;
		}
		position++;		/* 下一个位置 */
		if (i<0)
		{
			break;
		}
	}

	/* 第二行 */
	high_num = g_row2;
	position = 1;			/* 第几个8进制数 */
	while (high_num>0)
	{
		low_num = high_num%10;		/* 得到这个8进制数 */
		if (low_num>7)				/* 判断是否为8进制数 */
		{
			return 0;
		}
		high_num = high_num/10;		/* 记下其余部分 */

		/* 将8进制数转为二进制并保存 */
		for (i=k_column-(position-1)*3-1; i>=0 && i>=k_column-position*3; i--)
		{
			*(mx_g_turbo+k_column+i) = low_num%2;
			low_num = low_num/2;
		}
		position++;					/* 下一个位置 */
		if (i<0)
		{
			break;
		}
	}
	return 1;
}

/*---------------------------------------------------------------
函数:
	void int2bin(int intstat, int *tempstat, int length)
介绍:
	十进制数转为二进制序列.
参数:
	输入参数:
		intstat - 十进制数.
		length - 要得到的二进制序列长度．
	输出参数:
		bin_stat - 二进制序列首址.
返回值:
	无
---------------------------------------------------------------*/
void int2bin(int intstat, int *bin_stat, int length)
{
	int i, temp;

	temp = intstat;

	/* 除以2求余数 */
	for (i=length-1; i>=0; i--)
	{
		*(bin_stat+i) = temp%2;
		temp = temp/2;
	}
}

/*---------------------------------------------------------------
函数:
	int bin2int(int *binseq, int length)
介绍:
	二进制序列转为十进制数.
参数:
	输入参数:
		binseq - 二进制序列首址.
		length - 二进制序列长度．
	输出参数:
		无
返回值:
	得到的十进制数．
---------------------------------------------------------------*/
int bin2int(int *binseq, int length)
{
	int i, j, temp;
	int sum = 0;

	for (i=0; i<length; i++)
	{
		temp = 1;

		/* 计算该位权值 */
		for (j=1; j<=i; j++)
		{
			temp = temp * 2;
		}
		sum = sum + temp * (*(binseq+length-1-i));
	}

	return sum;
}

/*---------------------------------------------------------------
函数:
	int encode_bit(int inbit, int *stat)
介绍:
	比特编码器.
参数:
	输入参数:
		inbit -输入比特.
		stat - 当前寄存器状态首址.
	输出参数:
		stat - 编码后寄存器状态首址.
返回值:
	输出比特．
---------------------------------------------------------------*/
int encode_bit(int inbit, int *stat)
{
	int j;			/* 循环变量 */
	int output;		/* 输出比特 */

	/* 计算输出比特 */
	output = (*(turbo_g.g_matrix+turbo_g.K_num_col+0)) * inbit;

	for (j=1; j<turbo_g.K_num_col; j++)
	{
		output = (output + (*(turbo_g.g_matrix+turbo_g.K_num_col+j)) * (*(stat+j-1)))%2;
	}

	/* 修改状态序列 */
	for (j=turbo_g.K_num_col-2; j>0; j--)
	{
		*(stat+j)=*(stat+j-1);
	}

	*(stat+0) = inbit;

	return output;
}

/*---------------------------------------------------------------
函数:
	void gen_trellis()
介绍:
	生成Trellis.
参数:
	无
返回值:
	无
---------------------------------------------------------------*/
void gen_trellis()
{
	int i, j, k;						/* 循环变量 */
	int dk_turbo, ak_turbo, outbit;		/* 编码器内部比特和输出比特 */

	int *tempstat;						/* 状态序列 */

	if ((tempstat=(int *)malloc(sizeof(int)*M_num_reg))==NULL)
	{
	  printf("\n fail to allocate memory of tempstat \n");
	  exit(1);  
	}
	/* 生成后向输出和后向状态矩阵 */
	for (i=0; i<n_states; i++)			/* 状态循环 */
	{
		for (j=0; j<2; j++)				/* 输入为0,1 */
		{
			int2bin(i, tempstat, M_num_reg);	/* 将状态转为二进制序列 */

			/* dk:原始输入比特 */
			dk_turbo = j;

			/* 计算ak:叠加反馈后的输入比特，根据生成矩阵第一行（13） */
			ak_turbo = (*(turbo_g.g_matrix+0)) * dk_turbo;
			for (k=1; k<turbo_g.K_num_col; k++)
			{
				ak_turbo = ak_turbo + (*(turbo_g.g_matrix+k)) * (*(tempstat+k-1));
			}

			ak_turbo = ak_turbo % 2;

			/* 计算输出比特,修改状态序列，根据生成矩阵第二行（15） */
			outbit = encode_bit(ak_turbo, tempstat);

			/* 写入后向输出和后向状态矩阵 */
			*(turbo_trellis.mx_nextout+i*4+2*j)=2*dk_turbo-1;
			*(turbo_trellis.mx_nextout+i*4+2*j+1)=2*outbit-1;
			*(turbo_trellis.mx_nextstat+i*2+j)=bin2int(tempstat, M_num_reg);
		}			/* 输入循环结束 */

	}	/* 状态循环结束 */

	/* 生成前向输出和前向状态矩阵 */
	for (i=0; i<n_states; i++)	/* 状态循环 */
	{
		for (j=0; j<2; j++)/* 输入为0,1 */
		{
			*(turbo_trellis.mx_laststat+(*(turbo_trellis.mx_nextstat+i*2+j))*2+j) = i;
			*(turbo_trellis.mx_lastout+(*(turbo_trellis.mx_nextstat+i*2+j))*4+2*j)
				= *(turbo_trellis.mx_nextout+i*4+2*j);
			*(turbo_trellis.mx_lastout+(*(turbo_trellis.mx_nextstat+i*2+j))*4+2*j+1)
				= *(turbo_trellis.mx_nextout+i*4+2*j+1);
		}	/* 输入循环结束 */
	}	/* 状态循环结束 */

	free(tempstat);
}

/*---------------------------------------------------------------
函数:
	void TurboCodingInit()
介绍:
	Turbo码初始化函数.
参数:
	无
返回值:
	无
---------------------------------------------------------------*/
void TurboCodingInit()
{
	/* 初始化生成阵 */
	turbo_g.N_num_row = 3;				/* 行数 */
	turbo_g.K_num_col = COLUMN_OF_G;	/* 列数 */

	SYMBOL_NUM = (int)((double)length_after_code/(MODULATION));

	if ((turbo_g.g_matrix=(int *)malloc(turbo_g.N_num_row*turbo_g.K_num_col*sizeof(int)))==NULL)
	{
		printf("\n fail to allocate memory of turbo_g\n");
		exit(1);  
	}

	/* 得到生成阵 */
	if (!gen_g_matrix(turbo_g.K_num_col, G_ROW_1, G_ROW_2, turbo_g.g_matrix))
	{
		printf("error number of G\n");
		exit(1);
	}

	/* 生成Trellis */
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

	/* 为随机交织器分配内存 */
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

	/*生成QPP交织器下标*/
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
函数:
	void rsc_encode(int *source, int *rsc, int terminated, int len_info)
介绍:
	RSC编码器.
参数:
	输入参数:
		source -源数据序列首址.
		len_info - 源数据序列长度．
		terminated - 是否结尾: 1-结尾, 0-不结尾.
	输出参数:
		RSC - 编码后数据序列首址.
返回值:
	无．
---------------------------------------------------------------*/
void rsc_encode(int *source, int *rsc, int terminated, int len_info)
{
	int i, j;			/* 循环变量 */

	int *state;			/* 状态序列 */
	int dk_turbo, ak_turbo, outbit;		/* 编码器内的dk,ak和编码输出比特 */

	int len_total;						/* 总长度 */

	if ((state=(int *)malloc(M_num_reg*sizeof(int)))==NULL)
	{
	  printf("\n fail to allocate memory of state \n");
	  exit(1);  
	}

	/* 计算总长度 */
	len_total = len_info+M_num_reg;

	/* 初状态为0 */
	for (i=0; i<M_num_reg; i++)
	{
		*(state+i) = 0;
	}

	for (i=0; i<len_total; i++)							/* 逐比特编码 */
	{
		if (!terminated || (terminated && i<len_info))	/* 对信息比特 */
		{
			dk_turbo = *(source+i);
		}
		else											/* 结尾Trellis */
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

		/* 计算ak */
		ak_turbo = *(turbo_g.g_matrix+0) * dk_turbo;
		for (j=1; j<turbo_g.K_num_col; j++)
		{
			ak_turbo = ak_turbo + (*(turbo_g.g_matrix+j))*(*(state+j-1));
		}

		ak_turbo = ak_turbo%2;

		/* 对ak进行比特编码 */
		outbit = encode_bit(ak_turbo, state);

		/* 写dk和输出比特 */
		*(rsc+2*i) = dk_turbo;
		*(rsc+2*i+1) = outbit;
	}				/* 逐比特编码结束 */

	free(state);
}

/*---------------------------------------------------------------
函数:
	void encoderm_turbo(int *source, int *send_turbo, int len_info)

介绍:
	Turbo编码调制.
参数:
	输入参数:
		source -源数据序列首址.
		len_info - 源数据序列长度．
		type_flow - 信息类型: 1-业务信息, 0补充信息.
	输出参数:
		send_turbo - 编码调制后数据序列首址.
返回值:
	无．
---------------------------------------------------------------*/
void encoderm_turbo(int *source, int *send_turbo, int len_info)
{
	int i;									/* 循环变量 */
	int len_total = len_info + M_num_reg;	/* 总长度 */

	int *rsc1, *rsc2;		/* 两个RSC编码器的输出 */
	
	int *input2;			/* RSC2的输入 */

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

	/* 交织信源比特 */
	randominterleaver_int(source, input2, index_randomintlvr, len_info);

	/* RSC2 */
	rsc_encode(input2, rsc2, TERMINATED, len_info);

	/* 信息位复用调制 */
	for (i=0; i<len_info; i++)
	{
		*(send_turbo+3*i) = *(rsc1+2*i);
		*(send_turbo+3*i+1) = *(rsc1+2*i+1);
		*(send_turbo+3*i+2) = *(rsc2+2*i+1);
	}
	
	/* 结尾位复用调制 */
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
函数:
	double random_turbo()
介绍:
	生成0-1均匀分布随机数.
参数:
	无
返回值:
	生成的0-1均匀分布随机数．
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
		*(index+i) = (f1*i+(((f2*i)%source_length)*i)%source_length)%source_length; // 加类型转换是为了防止数值太大而溢出！
	}
}

/*---------------------------------------------------------------
函数:
	void gen_rand_index(int length, int type_flow)
介绍:
	生成随机交织器的下标.
参数:
	输入参数:
		length - 交织器长度.
		type_flow - 信息类型: 1-业务信息, 0-补充信息.
	输出参数:
		无
返回值:
	无
---------------------------------------------------------------*/
void gen_rand_index(int length, int *index)
{
	int *index_random;			/* 需先生成的0-1均匀分布随机数序列 */
	double *tempindex;			/* 用以选择随机交织首址 */
	double tempmax;				/* 最大值 */
	int selectedscr;			/* 选中的下标 */
	int i, j;					/* 循环变量 */

	if ((tempindex=(double *)malloc((length)*sizeof(double)))==NULL)
	{
	  printf("\n fail to allocate memory of tempindex \n");
	  exit(1);  
	}

	/* 由数据类型选随机交织器 */
	index_random = index;

	/* 生成的0-1均匀分布随机数序列 */
	for (i=0; i<length; i++)
	{
		*(tempindex+i) = random_turbo();
	}

	for (i=0; i<length; i++)	
	{
		/* 找到tempindex中的最大值对应的下标 */
		tempmax = 0.0;

		for (j=0; j<length; j++)
		{
			if (*(tempindex+j) >= tempmax)
			{
				tempmax = *(tempindex+j);
				selectedscr = j;
			}
		}

		/* 交织器该位置为该下标 */
		*(index_random+i) = selectedscr;

		/* tempindex该位赋0 */
		*(tempindex+selectedscr) = 0.0;
	}

	free(tempindex);
}

/*---------------------------------------------------------------
函数:
	TurboEnCoding(int *source, int *coded_source,
						int *source_length)
介绍:
	Turbo码业务信息编码函数.
参数:
	输入参数: source - 源bit序列首址.
			  source_length - 源bit序列长度.
	输出参数: coded_source - 编码后的序列首址.
返回值:
	无
---------------------------------------------------------------*/
void TurboEnCoding(int *source, int *coded_source, int source_length)
{
	int i;							/* 循环变量 */

	int *temp_send = NULL;			
	int *send = NULL;				/* 发送序列首址 */

	int length_info = source_length;		/* 信息位长度 */

	/* 申请内存 */
	if ((send=(int *)malloc((3*length_info+4*M_num_reg)*sizeof(int)))==NULL)
	{
		printf("\n fail to allocate memory of send \n");
		exit(1);
	}

	/* 生成随机交织器下标 */
//	gen_rand_index(length_info, index_randomintlvr);

	encoderm_turbo(source, send, length_info);	/* 编码 */

	temp_send = send;

	/* 写输出序列 */
	for (i=0; i<(3*length_info+4*M_num_reg); i++)
	{
		*(coded_source+i) = *(temp_send+i);
	}

	free(send);
}

/*---------------------------------------------------------------
函数:
	double get_max(double *data_seq, int length)
介绍:
	得到序列的最大值.
参数:
	输入参数:
		data_seq - 序列首址.
		length - 序列长度.
	输出参数:
		无
返回值:
	序列中的最大值．
---------------------------------------------------------------*/
double get_max(double *data_seq, int length)
{
	int i;		/* 循环变量 */
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
函数:
	double E_algorithm(double x, double y)
介绍:
	Log-MAP中的E算法:log(exp(x) + exp(y)) = max(x,y)+f(-|y-x|).
参数:
	输入参数:
		x,y - 输入数.
	输出参数:
		无
返回值:
	输出．
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
		/* 查找下标 */
		for (i=0; i<16 && temp>=lookup_index_Log_MAP[i]; i++)
		{
			;
		}
		/* 查找f(-|y-x|) */
		temp = (double)lookup_table_Log_MAP[i-1];
	}
	
	/* 返回max(x,y)+f(-|y-x|) */
	return ( (x>y?x:y) + temp );
}	

/*---------------------------------------------------------------
函数:
	double E_algorithm_seq(double *data_seq, int length)
介绍:
	序列的E算法.
参数:
	输入参数:
		data_seq - 序列首址.
		length - 序列长度
	输出参数:
		无
返回值:
	结果．
---------------------------------------------------------------*/
double E_algorithm_seq(double *data_seq, int length)
{
	int i;			/* 循环变量 */
	double temp;

	/* 每两个进行E算法,再与下一个进行E算法 */
	temp = E_algorithm(*(data_seq+0), *(data_seq+1));
	for (i=2; i<length; i++)
	{
		temp = E_algorithm(temp, *(data_seq+i));
	}
	return temp;
}


/*---------------------------------------------------------------
函数:
	void decision(double *LLR_seq, int length, int *output)
介绍:
	判决.
参数:
	输入参数:
		LLR_seq - LLR序列首址.
		length - 序列长度
	输出参数:
		output - 判决后输出序列首址
返回值:
	无．
---------------------------------------------------------------*/


/*---------------------------------------------------------------
函数:
	void decision(double *LLR_seq, int length, int *output)
介绍:
	判决.
参数:
	输入参数:
		LLR_seq - LLR序列首址.
		length - 序列长度
	输出参数:
		output - 判决后输出序列首址
返回值:
	无．
---------------------------------------------------------------*/
void decision(double *LLR_seq, int length, int *output)
{
	int i;

	for (i=0; i<length; i++)
	{
		/* 小于0判0 */
		if (*(LLR_seq+i) < 0)
		{
			*(output+i) = 0;
		}
		/* 大于0判1 */
		else
		{
			*(output+i) = 1;
		}
	}		
}

/*---------------------------------------------------------------
函数:
	void Log_MAP_decoder(double *recs_turbo, double *La_turbo, int terminated, 
						double *LLR_all_turbo, int len_total)
介绍:
	Log_MAP译码器.
参数:
	输入参数:
		rec_turbo -乘Lc后的接收序列首址.
		La_turbo - 外信息．
		terminated - 是否结尾.
		len_total - 输入序列长度.
	输出参数:
		LLR_all_turbo - 似然比.
返回值:
	无．
---------------------------------------------------------------*/
void Log_MAP_decoder(double *recs_turbo, double *La_turbo, int terminated, double *LLR_all_turbo, int len_total)
{
	int i, j;			/* 循环变量 */

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
	
	/*===== 初始化alpha和beta =====*/
	*(alpha_Log+0) = 0;
	*(beta_Log+len_total) = 0;

	for (i=1; i<n_states; i++)
	{
		*(alpha_Log+i*(len_total+1)+0) = (double)-INFTY;

		/* 结尾 */
		if (terminated)
		{
			*(beta_Log+i*(len_total+1)+len_total) = (double)-INFTY;
		}
		/* 不结尾 */
		else
		{
			*(beta_Log+i*(len_total+1)+len_total) = 0;
		}
	}
	/*======== 计算Gama ========*/
	for (i=0; i<len_total; i++)		/* 逐比特循环 */
	{
		for (j=0; j<n_states; j++)	/* 按状态循环 */
		{
			/* 见论文gama计算式 */
			*(gama_Log+j*len_total*2+i*2+0) 
				= -*(recs_turbo+2*i) + *(recs_turbo+2*i+1)*(*(turbo_trellis.mx_nextout+j*4+1)) - *(La_turbo+i)/2;
			*(gama_Log+j*len_total*2+i*2+1)
				= *(recs_turbo+2*i) + *(recs_turbo+2*i+1)*(*(turbo_trellis.mx_nextout+j*4+3)) + *(La_turbo+i)/2;
		}	/* 按状态循环结束 */
	}		/* 逐比特循环结束 */

	/*======== 前向计算alpha ========*/
	for (i=1; i<len_total+1; i++)	/* 逐比特循环 */
	{
		for (j=0; j<n_states; j++)	/* 按状态循环 */
		{
			tempx = *(gama_Log+(*(turbo_trellis.mx_laststat+j*2+0))*len_total*2+(i-1)*2+0)
					+ *(alpha_Log+(*(turbo_trellis.mx_laststat+j*2+0))*(len_total+1)+i-1);
			tempy = *(gama_Log+(*(turbo_trellis.mx_laststat+j*2+1))*len_total*2+(i-1)*2+1)
					+ *(alpha_Log+(*(turbo_trellis.mx_laststat+j*2+1))*(len_total+1)+i-1);
			*(alpha_Log+j*(len_total+1)+i) = E_algorithm(tempx, tempy);
		}	/* 按状态循环结束 */

		/* 计算tempmax[i],用于规一化alpha和beta */
		for (j=0; j<n_states; j++)
		{
			if (*(tempmax+i) < *(alpha_Log+j*(len_total+1)+i))
			{
				*(tempmax+i) = *(alpha_Log+j*(len_total+1)+i);
			}
		}

		/* 规一化alpha */
		for(j=0; j<n_states; j++)
		{
			*(alpha_Log+j*(len_total+1)+i) = *(alpha_Log+j*(len_total+1)+i) - *(tempmax+i);
		}

	}	/* 逐比特循环结束 */

	/*======== 反向计算beta ========*/
	for (i=len_total-1; i>=0; i--)	/* 逐比特循环 */
	{
		for (j=0; j<n_states; j++)	/* 按状态循环 */
		{
			tempx = *(gama_Log+j*len_total*2+i*2+0) 
					+ *(beta_Log+(*(turbo_trellis.mx_nextstat+j*2+0))*(len_total+1)+i+1);
			tempy = *(gama_Log+j*len_total*2+i*2+1) 
					+ *(beta_Log+(*(turbo_trellis.mx_nextstat+j*2+1))*(len_total+1)+i+1);

			*(beta_Log+j*(len_total+1)+i) = E_algorithm(tempx, tempy);
		}	/* 按状态循环结束 */

		/* 规一化beta */
		for (j=0; j<n_states; j++)
		{
			*(beta_Log+j*(len_total+1)+i) = *(beta_Log+j*(len_total+1)+i) - *(tempmax+i+1);
		}
	}	/* 逐比特循环结束 */

	/*=== 计算似然比LLR ===*/
	for (i=0; i<len_total; i++)		/* 逐比特循环 */
	{
		for (j=0; j<n_states; j++)	/* 按状态循环 */
		{
			*(temp0+j) = *(gama_Log+(*(turbo_trellis.mx_laststat+j*2+0))*len_total*2+i*2+0) 
				+ *(alpha_Log+*(turbo_trellis.mx_laststat+j*2+0)*(len_total+1)+i)
				+ *(beta_Log+ j*(len_total+1)+i+1);

			*(temp1+j) = *(gama_Log+(*(turbo_trellis.mx_laststat+j*2+1))*len_total*2+i*2+1) 
				+ *(alpha_Log+*(turbo_trellis.mx_laststat+j*2+1)*(len_total+1)+i)
				+ *(beta_Log+j*(len_total+1)+i+1);
		}	/* 按状态循环结束 */

		/* 计算似然比 */
		*(LLR_all_turbo+i) = E_algorithm_seq(temp1, n_states) - E_algorithm_seq(temp0, n_states);
	}	/* 逐比特循环结束 */

	free(alpha_Log);
	free(beta_Log);
	free(gama_Log);
	free(tempmax);
	free(temp0);
	free(temp1);
}

/*---------------------------------------------------------------
函数:
	void Log_MAP_decoder_fix(int *recs_turbo, long *La_turbo, int terminated, 
						long *LLR_all_turbo, int len_total)
介绍:
	Log_MAP译码器.
参数:
	输入参数:
		rec_turbo -乘Lc后的接收序列首址.
		La_turbo - 外信息．
		terminated - 是否结尾.
		len_total - 输入序列长度.
	输出参数:
		LLR_all_turbo - 似然比.
返回值:
	无．
---------------------------------------------------------------*/



/*---------------------------------------------------------------
函数:
	void demultiplex(double *rec_turbo, int len_info, double *yk_turbo)
介绍:
	解复用.
参数:
	输入参数:
		rec_turbo -接收数据序列首址.
		len_info - 接收数据序列长度．
	输出参数:
		yk_turbo - 解复用后数据序列首址.
返回值:
	无．
---------------------------------------------------------------*/
void demultiplex(double *rec_turbo, int len_info, double *yk_turbo)
{
	int i;			/* 循环变量 */

	int len_total = len_info+M_num_reg;		/* 总长度 */

	double *info2, *inted_info2;			/* 信息位和交织后的信息位 */

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
	
	/* 对信息比特 */
	for(i=0; i<len_info; i++)
	{
		*(info2+i) = *(yk_turbo+2*i) = *(rec_turbo+3*i);
		*(yk_turbo+2*i+1) = *(rec_turbo+3*i+1);
		*(yk_turbo+2*len_total+2*i+1) = *(rec_turbo+3*i+2);
	}

	/* 交织信息比特 */
	randominterleaver_double(info2, inted_info2, index_randomintlvr, len_info);

	for (i=0; i<len_info; i++)
	{
		*(yk_turbo+2*len_total+2*i) = *(inted_info2+i);
	}

	/* 对尾比特 */
	for (i=0; i<2*M_num_reg; i++)
	{
		*(yk_turbo+2*len_info+i) = *(rec_turbo+3*len_info+i);
		*(yk_turbo+2*len_total+2*len_info+i) = *(rec_turbo+3*len_info+2*M_num_reg+i);
	}

	free(info2);
	free(inted_info2);
}




/*---------------------------------------------------------------
函数:
	void TurboDecoding(double *flow_for_decode, int *flow_decoded,
						  int *flow_length, double EbN0dB)
介绍:
	Turbo码业务信息译码函数.
参数:
	输入参数: flow_for_decode - 输入序列首址.
			  flow_length - 输入序列长度.
			  EbN0dB - 信噪比.
	输出参数: flow_decoded - 译码后的序列首址.
返回值:
	无
---------------------------------------------------------------*/
void TurboDecoding(double *flow_for_decode, int *flow_decoded,int flow_length)
{
	int i;							/* 循环变量 */
	int length_info, length_total;	/* 信息位长和总长度 */
	int iteration;					/* 叠代循环变量 */

	float *receive_punc = NULL;					/* 接收数据 */ 
	double *yk_turbo = NULL;					/* 包含ys和yp */
	double *La_turbo, *Le_turbo, *LLR_all_turbo;		/* 外信息和似然比 */

	int *tempout;


	/* 计算信息位长度 */
	length_info = (flow_length-4*M_num_reg)/3;

	/* 总长度 */
	length_total = length_info+M_num_reg;

	/* 申请内存 */
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

	/* 解复用 */
	demultiplex(flow_for_decode, length_info, yk_turbo);	/* 解复用 */
	
	/* 初始化外信息和似然信息 */
	for (i=0; i<length_total; i++)
	{
		*(La_turbo+i) = *(Le_turbo+i) = *(LLR_all_turbo+i) = 0;
	}

	for (iteration=0; iteration<N_ITERATION; iteration++)		/* 开始叠代 */
	{
		/* 译码器1: */
		/* 解交织来自译码器2的外信息 */
		random_deinterlvr_double(La_turbo, Le_turbo, index_randomintlvr, length_info);
		
		/* 结尾处赋0 */
		for (i=length_info; i<length_total; i++)
		{
			*(La_turbo+i) = 0;
		}

			/* Log-MAP算法 */
			Log_MAP_decoder(yk_turbo, La_turbo, TERMINATED, LLR_all_turbo, length_total);
		

		/* 计算外信息 */
		for (i=0; i<length_total; i++)
		{
			//csqu-070717
			*(Le_turbo+i) = *(LLR_all_turbo+i) - *(La_turbo+i) - 2*(*(yk_turbo+2*i));
		}

		/* 译码器2: */
		/* 交织来自译码器1的外信息 */
		randominterleaver_double(Le_turbo, La_turbo, index_randomintlvr, length_info);

		/* 结尾处赋0 */
		for (i=length_info; i<length_total; i++)
		{
			*(La_turbo+i) = 0;
		}
						
			/* Log-MAP算法 */
			Log_MAP_decoder(yk_turbo+2*length_total, La_turbo, TERMINATED, LLR_all_turbo, length_total);
		

		/* 计算外信息 */
		for(i=0; i<length_total; i++)
		{
			//csqu-070717
			*(Le_turbo+i) = *(LLR_all_turbo+i) - *(La_turbo+i) - 2*(*(yk_turbo+2*length_total+2*i));
		}
		/* 译码器2结束 */
		decision(LLR_all_turbo, length_total, tempout);

	/* 对判决比特解交织 */
		random_deinterlvr_int(flow_decoded+length_info*iteration, tempout, index_randomintlvr, length_info);
	}
	
	/* 对全部似然比判决 */

	


	free(receive_punc);
	free(yk_turbo);

	free(La_turbo);
	free(Le_turbo);
	free(LLR_all_turbo);

	free(tempout);
}


void dectobin(double *flow_for_change, int *flow_changed, int flow_len, int integer_len, int decimal_len)
{
    int i;/* 循环变量 */
	int total_len = integer_len + decimal_len;

		/* 申请内存 */
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
	flow_changed_f<< "-----flow_changed文件开始-----" << endl;
	for (i=0; i<flow_len; i++) flow_changed_f << flow_changed[i] << endl;

	flow_changed_f<< "-----flow_changed文件结束-----" << endl;*/


}



void TurboCodingRelease()
{
	/* 释放生成阵内存 */
	free(turbo_g.g_matrix);

	/* 释放trellis结构内存 */
	free(turbo_trellis.mx_lastout);
	free(turbo_trellis.mx_laststat);
	free(turbo_trellis.mx_nextout);
	free(turbo_trellis.mx_nextstat);

	/* 释放随机交织器内存 */
	free(index_randomintlvr);
	free(index_randomintlvr2);

}
/*---------------------------------------------------------------
函数:
	void mgrns(double mean,double sigma,double seed,int n,double *a)
介绍:
	产生长度为n的高斯随机序列.
参数:
	输入参数:	mean - 均值
				sigma - 标准差
				seed - 一个随机种子
	输出参数:	a - 长度为n的高斯随机序列.
返回值:
	无
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
        }/*按照中心极限定理产生服从高斯分布的随机数*/
        *(a+k)=mean+(double)(sigma*(t-6.0));
    }
    return;
}
/*---------------------------------------------------------------
函数:
	void AWGN(int *send, double *r, double sigma, int totallength)
介绍:
	AWGN信道.
参数:
	输入参数: send - int型发送数据首址.
			  sigma - 噪声标准差.
			  totallength - 序列长度.
	输出参数: r - 经AWGN信道后的数据序列的首址.
返回值:
	无
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


