#include "main.h"
#include <fstream>
#include <iostream>
using namespace std;

extern int M_num_reg;
int main()
{
/*void main(int argc, char *argv[])
{
          int RB=atoi(argv[1]);
	int MCS=atoi(argv[2]);
	MODULATION = atoi(argv[3]);	//6;			//调制阶数：1,2,3,4,6
	length_after_code = atoi(argv[4]);//SYMBOL_NUM*MODULATION;	//编码凿孔后的长度	
	SYMBOL_NUM=length_after_code/MODULATION;//  atoi(argv[4]);//672

	source_length =atoi(argv[5]);//2688;			//输入信源长度
	f1 =atoi(argv[6]);// 127;
	f2 =atoi(argv[7]);//504; 					//与码长有关的交织器参数

	double EbN0start = atof(argv[8]);//0;		//仿真起始信噪比
	double EbN0end = atof(argv[9]);//4;		//最大仿真终止信噪比
	double EbN0step =atof(argv[10]);;		//仿真信噪比步长
 */
    
	
	//int RB=1;
	//int MCS=atoi(argv[2]);
	MODULATION = 1;			//调制阶数：1,2,3,4,6
	source_length = 6144;			//输入信源长度
	length_after_code = 6144*3+12;//SYMBOL_NUM*MODULATION;	//编码凿孔后的长度	
	//SYMBOL_NUM=length_after_code/MODULATION;//  atoi(argv[4]);//672
	SYMBOL_NUM=length_after_code;

	
	f1 =263;
	f2 =480; 					//与码长有关的交织器参数

	double EbN0start = 0;//0;		//仿真起始信噪比
	double EbN0end = 1;//4;		//最大仿真终止信噪比
	double EbN0step =0.1;		//仿真信噪比步长
	int FrameNum=100000;	

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

	int N_ITERATION=15;
	int temp[15],err_bit_num[15], err_block_num[15];
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
			
				if(err_block_num[N_ITERATION-1]>=50)//错够1000个块，跳出
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
	//		break;//误块率小于0.01，跳出
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
