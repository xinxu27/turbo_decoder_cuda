#include<iostream>
#include<stdio.h>
#include<stdlib.h>

//=============================================================================================================//

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
