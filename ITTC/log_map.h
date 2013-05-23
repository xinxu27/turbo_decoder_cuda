/*----------------------------------------------------------
* Copyright (c) 2003, 北京邮电大学移动通信实验室
* All rights reserved.
*
* 文件名称：turbo_code_Log_MAP.h
* 文件标识：
* 摘    要：Turbo编译码头文件.
*
* 当前版本：1.0
* 作    者：张鹏
* 完成日期：2003年12月1日
----------------------------------------------------------*/

#ifndef	TRUBO_CODE_LOG_MAP_H
#define	TRUBO_CODE_LOG_MAP_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

/*==================================================*/
/* 仿真参数 */
#define TERMINATED          1               /* 0不结尾；1结尾 */

#define TYPE_DECODER		1				/* 译码器类型:1-LogMAP译码
														  2-MAX-LogMAP译码
														  3-SOVA译码
                                                          4-const_LogMAP */
#define N_ITERATION			8				/* 译码叠代次数 */
#define MAX_FRAME_LENGTH	10000			/* 最大帧长 */
/*==================================================*/
/* 生成阵参数 */
#define COLUMN_OF_G		4					/* 生成阵列数 */
#define G_ROW_1			13                  /* 反馈抽头 */					 
#define G_ROW_2			15					/* 输出抽头 */

/*==================================================*/
/*==================================================*/

/* 定义数据结构 */

/* 生成阵结构 */
typedef struct
{
	int N_num_row;							/* 生成阵行数 */
	int K_num_col;							/* 生成阵列数 */
	int *g_matrix;							/* 生成阵首址 */
} TURBO_G;

/* Trellis结构 */
/*
mx_nextout(mx_lastout)为后(前)向输出数组首址:
行数:状态数, 列数:4
每行第一列和第三列为后(前)向的输入(1或-1),第二列和第四列为与之对应的输出(1或-1).

mx_nextstat(mx_laststat)为后(前)向状态数组首址:
行数:状态数, 列数:2
各列表示输入为1(0)时对应的后(前)向状态.
*/
typedef struct
{
	int *mx_nextout;		/* 后向输出矩阵 */						
	int *mx_nextstat;		/* 后向状态矩阵 */
	int *mx_lastout;		/* 前向输出矩阵 */
	int *mx_laststat;		/* 前向状态矩阵 */
	
} TURBO_TRELLIS;

/*==================================================*/
/* 无穷大 */
#ifndef INFTY
#define INFTY 1E20
#endif
/*==================================================*/

/* 全局变量 */

int *index_randomintlvr;		/* 业务信息随机交织器下标 */

int *index_randomintlvr2;

TURBO_G turbo_g;				/* 生成阵 */

TURBO_TRELLIS turbo_trellis;	/* Tellis结构 */	

double rate_coding;				/* 编码速率 */

void gen_qpp_index(int length, int *index);


#endif