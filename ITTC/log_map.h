/*----------------------------------------------------------
* Copyright (c) 2003, �����ʵ��ѧ�ƶ�ͨ��ʵ����
* All rights reserved.
*
* �ļ����ƣ�turbo_code_Log_MAP.h
* �ļ���ʶ��
* ժ    Ҫ��Turbo������ͷ�ļ�.
*
* ��ǰ�汾��1.0
* ��    �ߣ�����
* ������ڣ�2003��12��1��
----------------------------------------------------------*/

#ifndef	TRUBO_CODE_LOG_MAP_H
#define	TRUBO_CODE_LOG_MAP_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

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


#endif