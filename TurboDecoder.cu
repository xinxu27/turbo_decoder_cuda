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

#define MIN 1E-300
#define L_TOTAL 4096// if u want to use block interleave,L_TOTAL must = x^2
#define MAXITER 10
#define	FRAME_NUM 10
#define AlphaBetaTHREAD_NUM 4

#define THREAD_NUM 1024
#define BLOCK_NUM 4

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




#define M	3	// register length,=tail length
#define NSTATE	8	// = M^2
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
{	0,0,1,1,1,1,0,0,
	0,0,1,1,1,1,0,0
};
static const char NextOut[2][NSTATE] = // check bit based on current and input bit
{	-1,-1,1,1,1,1,-1,-1,
	-1,-1,1,1,1,1,-1,-1
};
// NextState[bk][current state]
static const BYTE NextState[2][NSTATE] = // next state based on current and input bit
{	0,4,5,1,2,6,7,3,
	4,0,1,5,6,2,3,7
};
// LastOut[bk][current state]
static const char LastOut[2][NSTATE] =	// trellis last check bit
{	-1,1,1,-1,-1,1,1,-1,
	-1,1,1,-1,-1,1,1,-1
};
// LastState[bk][current state]
static const BYTE LastState[2][NSTATE] =	// last state lead to current state by input bk
{	0,3,4,7,1,2,5,6,
	1,2,5,6,0,3,4,7
};
// TailBit[current state]
static const char TailBit[NSTATE] = // tail info bits when trellis is terminating
{	0,1,1,0,0,1,1,0
};



UINT m_Inter_table[L_TOTAL] = 
{0,95,318,669,1148,1755,2490,3353,248,1367,
2614,3989,1396,3027,690,2577,496,2639,814,3213,
1644,203,2986,1801,744,3911,3110,2437,1892,1475,
1186,1025,992,1087,1310,1661,2140,2747,3482,249,
1240,2359,3606,885,2388,4019,1682,3569,1488,3631,
1806,109,2636,1195,3978,2793,1736,807,6,3429,
2884,2467,2178,2017,1984,2079,2302,2653,3132,3739,
378,1241,2232,3351,502,1877,3380,915,2674,465,
2480,527,2798,1101,3628,2187,874,3785,2728,1799,
998,325,3876,3459,3170,3009,2976,3071,3294,3645,
28,635,1370,2233,3224,247,1494,2869,276,1907,
3666,1457,3472,1519,3790,2093,524,3179,1866,681,
3720,2791,1990,1317,772,355,66,4001,3968,4063,
190,541,1020,1627,2362,3225,120,1239,2486,3861,
1268,2899,562,2449,368,2511,686,3085,1516,75,
2858,1673,616,3783,2982,2309,1764,1347,1058,897,
864,959,1182,1533,2012,2619,3354,121,1112,2231,
3478,757,2260,3891,1554,3441,1360,3503,1678,4077,
2508,1067,3850,2665,1608,679,3974,3301,2756,2339,
2050,1889,1856,1951,2174,2525,3004,3611,250,1113,
2104,3223,374,1749,3252,787,2546,337,2352,399,
2670,973,3500,2059,746,3657,2600,1671,870,197,
3748,3331,3042,2881,2848,2943,3166,3517,3996,507,
1242,2105,3096,119,1366,2741,148,1779,3538,1329,
3344,1391,3662,1965,396,3051,1738,553,3592,2663,
1862,1189,644,227,4034,3873,3840,3935,62,413,
892,1499,2234,3097,4088,1111,2358,3733,1140,2771,
434,2321,240,2383,558,2957,1388,4043,2730,1545,
488,3655,2854,2181,1636,1219,930,769,736,831,
1054,1405,1884,2491,3226,4089,984,2103,3350,629,
2132,3763,1426,3313,1232,3375,1550,3949,2380,939,
3722,2537,1480,551,3846,3173,2628,2211,1922,1761,
1728,1823,2046,2397,2876,3483,122,985,1976,3095,
246,1621,3124,659,2418,209,2224,271,2542,845,
3372,1931,618,3529,2472,1543,742,69,3620,3203,
2914,2753,2720,2815,3038,3389,3868,379,1114,1977,
2968,4087,1238,2613,20,1651,3410,1201,3216,1263,
3534,1837,268,2923,1610,425,3464,2535,1734,1061,
516,99,3906,3745,3712,3807,4030,285,764,1371,
2106,2969,3960,983,2230,3605,1012,2643,306,2193,
112,2255,430,2829,1260,3915,2602,1417,360,3527,
2726,2053,1508,1091,802,641,608,703,926,1277,
1756,2363,3098,3961,856,1975,3222,501,2004,3635,
1298,3185,1104,3247,1422,3821,2252,811,3594,2409,
1352,423,3718,3045,2500,2083,1794,1633,1600,1695,
1918,2269,2748,3355,4090,857,1848,2967,118,1493,
2996,531,2290,81,2096,143,2414,717,3244,1803,
490,3401,2344,1415,614,4037,3492,3075,2786,2625,
2592,2687,2910,3261,3740,251,986,1849,2840,3959,
1110,2485,3988,1523,3282,1073,3088,1135,3406,1709,
140,2795,1482,297,3336,2407,1606,933,388,4067,
3778,3617,3584,3679,3902,157,636,1243,1978,2841,
3832,855,2102,3477,884,2515,178,2065,4080,2127,
302,2701,1132,3787,2474,1289,232,3399,2598,1925,
1380,963,674,513,480,575,798,1149,1628,2235,
2970,3833,728,1847,3094,373,1876,3507,1170,3057,
976,3119,1294,3693,2124,683,3466,2281,1224,295,
3590,2917,2372,1955,1666,1505,1472,1567,1790,2141,
2620,3227,3962,729,1720,2839,4086,1365,2868,403,
2162,4049,1968,15,2286,589,3116,1675,362,3273,
2216,1287,486,3909,3364,2947,2658,2497,2464,2559,
2782,3133,3612,123,858,1721,2712,3831,982,2357,
3860,1395,3154,945,2960,1007,3278,1581,12,2667,
1354,169,3208,2279,1478,805,260,3939,3650,3489,
3456,3551,3774,29,508,1115,1850,2713,3704,727,
1974,3349,756,2387,50,1937,3952,1999,174,2573,
1004,3659,2346,1161,104,3271,2470,1797,1252,835,
546,385,352,447,670,1021,1500,2107,2842,3705,
600,1719,2966,245,1748,3379,1042,2929,848,2991,
1166,3565,1996,555,3338,2153,1096,167,3462,2789,
2244,1827,1538,1377,1344,1439,1662,2013,2492,3099,
3834,601,1592,2711,3958,1237,2740,275,2034,3921,
1840,3983,2158,461,2988,1547,234,3145,2088,1159,
358,3781,3236,2819,2530,2369,2336,2431,2654,3005,
3484,4091,730,1593,2584,3703,854,2229,3732,1267,
3026,817,2832,879,3150,1453,3980,2539,1226,41,
3080,2151,1350,677,132,3811,3522,3361,3328,3423,
3646,3997,380,987,1722,2585,3576,599,1846,3221,
628,2259,4018,1809,3824,1871,46,2445,876,3531,
2218,1033,4072,3143,2342,1669,1124,707,418,257,
224,319,542,893,1372,1979,2714,3577,472,1591,
2838,117,1620,3251,914,2801,720,2863,1038,3437,
1868,427,3210,2025,968,39,3334,2661,2116,1699,
1410,1249,1216,1311,1534,1885,2364,2971,3706,473,
1464,2583,3830,1109,2612,147,1906,3793,1712,3855,
2030,333,2860,1419,106,3017,1960,1031,230,3653,
3108,2691,2402,2241,2208,2303,2526,2877,3356,3963,
602,1465,2456,3575,726,2101,3604,1139,2898,689,
2704,751,3022,1325,3852,2411,1098,4009,2952,2023,
1222,549,4,3683,3394,3233,3200,3295,3518,3869,
252,859,1594,2457,3448,471,1718,3093,500,2131,
3890,1681,3696,1743,4014,2317,748,3403,2090,905,
3944,3015,2214,1541,996,579,290,129,96,191,
414,765,1244,1851,2586,3449,344,1463,2710,4085,
1492,3123,786,2673,592,2735,910,3309,1740,299,
3082,1897,840,4007,3206,2533,1988,1571,1282,1121,
1088,1183,1406,1757,2236,2843,3578,345,1336,2455,
3702,981,2484,19,1778,3665,1584,3727,1902,205,
2732,1291,4074,2889,1832,903,102,3525,2980,2563,
2274,2113,2080,2175,2398,2749,3228,3835,474,1337,
2328,3447,598,1973,3476,1011,2770,561,2576,623,
2894,1197,3724,2283,970,3881,2824,1895,1094,421,
3972,3555,3266,3105,3072,3167,3390,3741,124,731,
1466,2329,3320,343,1590,2965,372,2003,3762,1553,
3568,1615,3886,2189,620,3275,1962,777,3816,2887,
2086,1413,868,451,162,1,4064,63,286,637,
1116,1723,2458,3321,216,1335,2582,3957,1364,2995,
658,2545,464,2607,782,3181,1612,171,2954,1769,
712,3879,3078,2405,1860,1443,1154,993,960,1055,
1278,1629,2108,2715,3450,217,1208,2327,3574,853,
2356,3987,1650,3537,1456,3599,1774,77,2604,1163,
3946,2761,1704,775,4070,3397,2852,2435,2146,1985,
1952,2047,2270,2621,3100,3707,346,1209,2200,3319,
470,1845,3348,883,2642,433,2448,495,2766,1069,
3596,2155,842,3753,2696,1767,966,293,3844,3427,
3138,2977,2944,3039,3262,3613,4092,603,1338,2201,
3192,215,1462,2837,244,1875,3634,1425,3440,1487,
3758,2061,492,3147,1834,649,3688,2759,1958,1285,
740,323,34,3969,3936,4031,158,509,988,1595,
2330,3193,88,1207,2454,3829,1236,2867,530,2417,
336,2479,654,3053,1484,43,2826,1641,584,3751,
2950,2277,1732,1315,1026,865,832,927,1150,1501,
1980,2587,3322,89,1080,2199,3446,725,2228,3859,
1522,3409,1328,3471,1646,4045,2476,1035,3818,2633,
1576,647,3942,3269,2724,2307,2018,1857,1824,1919,
2142,2493,2972,3579,218,1081,2072,3191,342,1717,
3220,755,2514,305,2320,367,2638,941,3468,2027,
714,3625,2568,1639,838,165,3716,3299,3010,2849,
2816,2911,3134,3485,3964,475,1210,2073,3064,87,
1334,2709,116,1747,3506,1297,3312,1359,3630,1933,
364,3019,1706,521,3560,2631,1830,1157,612,195,
4002,3841,3808,3903,30,381,860,1467,2202,3065,
4056,1079,2326,3701,1108,2739,402,2289,208,2351,
526,2925,1356,4011,2698,1513,456,3623,2822,2149,
1604,1187,898,737,704,799,1022,1373,1852,2459,
3194,4057,952,2071,3318,597,2100,3731,1394,3281,
1200,3343,1518,3917,2348,907,3690,2505,1448,519,
3814,3141,2596,2179,1890,1729,1696,1791,2014,2365,
2844,3451,90,953,1944,3063,214,1589,3092,627,
2386,177,2192,239,2510,813,3340,1899,586,3497,
2440,1511,710,37,3588,3171,2882,2721,2688,2783,
3006,3357,3836,347,1082,1945,2936,4055,1206,2581,
4084,1619,3378,1169,3184,1231,3502,1805,236,2891,
1578,393,3432,2503,1702,1029,484,67,3874,3713,
3680,3775,3998,253,732,1339,2074,2937,3928,951,
2198,3573,980,2611,274,2161,80,2223,398,2797,
1228,3883,2570,1385,328,3495,2694,2021,1476,1059,
770,609,576,671,894,1245,1724,2331,3066,3929,
824,1943,3190,469,1972,3603,1266,3153,1072,3215,
1390,3789,2220,779,3562,2377,1320,391,3686,3013,
2468,2051,1762,1601,1568,1663,1886,2237,2716,3323,
4058,825,1816,2935,86,1461,2964,499,2258,49,
2064,111,2382,685,3212,1771,458,3369,2312,1383,
582,4005,3460,3043,2754,2593,2560,2655,2878,3229,
3708,219,954,1817,2808,3927,1078,2453,3956,1491,
3250,1041,3056,1103,3374,1677,108,2763,1450,265,
3304,2375,1574,901,356,4035,3746,3585,3552,3647,
3870,125,604,1211,1946,2809,3800,823,2070,3445,
852,2483,146,2033,4048,2095,270,2669,1100,3755,
2442,1257,200,3367,2566,1893,1348,931,642,481,
448,543,766,1117,1596,2203,2938,3801,696,1815,
3062,341,1844,3475,1138,3025,944,3087,1262,3661,
2092,651,3434,2249,1192,263,3558,2885,2340,1923,
1634,1473,1440,1535,1758,2109,2588,3195,3930,697,
1688,2807,4054,1333,2836,371,2130,4017,1936,4079,
2254,557,3084,1643,330,3241,2184,1255,454,3877,
3332,2915,2626,2465,2432,2527,2750,3101,3580,91,
826,1689,2680,3799,950,2325,3828,1363,3122,913,
2928,975,3246,1549,4076,2635,1322,137,3176,2247,
1446,773,228,3907,3618,3457,3424,3519,3742,4093,
476,1083,1818,2681,3672,695,1942,3317,724,2355,
18,1905,3920,1967,142,2541,972,3627,2314,1129,
72,3239,2438,1765,1220,803,514,353,320,415,
638,989,1468,2075,2810,3673,568,1687,2934,213,
1716,3347,1010,2897,816,2959,1134,3533,1964,523,
3306,2121,1064,135,3430,2757,2212,1795,1506,1345,
1312,1407,1630,1981,2460,3067,3802,569,1560,2679,
3926,1205,2708,243,2002,3889,1808,3951,2126,429,
2956,1515,202,3113,2056,1127,326,3749,3204,2787,
2498,2337,2304,2399,2622,2973,3452,4059,698,1561,
2552,3671,822,2197,3700,1235,2994,785,2800,847,
3118,1421,3948,2507,1194,9,3048,2119,1318,645,
100,3779,3490,3329,3296,3391,3614,3965,348,955,
1690,2553,3544,567,1814,3189,596,2227,3986,1777,
3792,1839,14,2413,844,3499,2186,1001,4040,3111,
2310,1637,1092,675,386,225,192,287,510,861,
1340,1947,2682,3545,440,1559,2806,85,1588,3219,
882,2769,688,2831,1006,3405,1836,395,3178,1993,
936,7,3302,2629,2084,1667,1378,1217,1184,1279,
1502,1853,2332,2939,3674,441,1432,2551,3798,1077,
2580,115,1874,3761,1680,3823,1998,301,2828,1387,
74,2985,1928,999,198,3621,3076,2659,2370,2209,
2176,2271,2494,2845,3324,3931,570,1433,2424,3543,
694,2069,3572,1107,2866,657,2672,719,2990,1293,
3820,2379,1066,3977,2920,1991,1190,517,4068,3651,
3362,3201,3168,3263,3486,3837,220,827,1562,2425,
3416,439,1686,3061,468,2099,3858,1649,3664,1711,
3982,2285,716,3371,2058,873,3912,2983,2182,1509,
964,547,258,97,64,159,382,733,1212,1819,
2554,3417,312,1431,2678,4053,1460,3091,754,2641,
560,2703,878,3277,1708,267,3050,1865,808,3975,
3174,2501,1956,1539,1250,1089,1056,1151,1374,1725,
2204,2811,3546,313,1304,2423,3670,949,2452,4083,
1746,3633,1552,3695,1870,173,2700,1259,4042,2857,
1800,871,70,3493,2948,2531,2242,2081,2048,2143,
2366,2717,3196,3803,442,1305,2296,3415,566,1941,
3444,979,2738,529,2544,591,2862,1165,3692,2251,
938,3849,2792,1863,1062,389,3940,3523,3234,3073,
3040,3135,3358,3709,92,699,1434,2297,3288,311,
1558,2933,340,1971,3730,1521,3536,1583,3854,2157,
588,3243,1930,745,3784,2855,2054,1381,836,419,
130,4065,4032,31,254,605,1084,1691,2426,3289,
184,1303,2550,3925,1332,2963,626,2513,432,2575,
750,3149,1580,139,2922,1737,680,3847,3046,2373,
1828,1411,1122,961,928,1023,1246,1597,2076,2683,
3418,185,1176,2295,3542,821,2324,3955,1618,3505,
1424,3567,1742,45,2572,1131,3914,2729,1672,743,
4038,3365,2820,2403,2114,1953,1920,2015,2238,2589,
3068,3675,314,1177,2168,3287,438,1813,3316,851,
2610,401,2416,463,2734,1037,3564,2123,810,3721,
2664,1735,934,261,3812,3395,3106,2945,2912,3007,
3230,3581,4060,571,1306,2169,3160,183,1430,2805,
212,1843,3602,1393,3408,1455,3726,2029,460,3115,
1802,617,3656,2727,1926,1253,708,291,2,3937,
3904,3999,126,477,956,1563,2298,3161,56,1175,
2422,3797,1204,2835,498,2385,304,2447,622,3021,
1452,11,2794,1609,552,3719,2918,2245,1700,1283,
994,833,800,895,1118,1469,1948,2555,3290,57,
1048,2167,3414,693,2196,3827,1490,3377,1296,3439,
1614,4013,2444,1003,3786,2601,1544,615,3910,3237,
2692,2275,1986,1825,1792,1887,2110,2461,2940,3547,
186,1049,2040,3159,310,1685,3188,723,2482,273,
2288,335,2606,909,3436,1995,682,3593,2536,1607,
806,133,3684,3267,2978,2817,2784,2879,3102,3453,
3932,443,1178,2041,3032,55,1302,2677,84,1715,
3474,1265,3280,1327,3598,1901,332,2987,1674,489,
3528,2599,1798,1125,580,163,3970,3809,3776,3871,
4094,349,828,1435,2170,3033,4024,1047,2294,3669,
1076,2707,370,2257,176,2319,494,2893,1324,3979,
2666,1481,424,3591,2790,2117,1572,1155,866,705,
672,767,990,1341,1820,2427,3162,4025,920,2039,
3286,565,2068,3699,1362,3249,1168,3311,1486,3885,
2316,875,3658,2473,1416,487,3782,3109,2564,2147,
1858,1697,1664,1759,1982,2333,2812,3419,58,921,
1912,3031,182,1557,3060,595,2354,145,2160,207,
2478,781,3308,1867,554,3465,2408,1479,678,5,
3556,3139,2850,2689,2656,2751,2974,3325,3804,315,
1050,1913,2904,4023,1174,2549,4052,1587,3346,1137,
3152,1199,3470,1773,204,2859,1546,361,3400,2471,
1670,997,452,35,3842,3681,3648,3743,3966,221,
700,1307,2042,2905,3896,919,2166,3541,948,2579,
242,2129,48,2191,366,2765,1196,3851,2538,1353,
296,3463,2662,1989,1444,1027,738,577,544,639,
862,1213,1692,2299,3034,3897,792,1911,3158,437,
1940,3571,1234,3121,1040,3183,1358,3757,2188,747,
3530,2345,1288,359,3654,2981,2436,2019,1730,1569,
1536,1631,1854,2205,2684,3291,4026,793,1784,2903,
54,1429,2932,467,2226,17,2032,79,2350,653,
3180,1739,426,3337,2280,1351,550,3973,3428,3011,
2722,2561,2528,2623,2846,3197,3676,187,922,1785,
2776,3895,1046,2421,3924,1459,3218,1009,3024,1071,
3342,1645,76,2731,1418,233,3272,2343,1542,869,
324,4003,3714,3553,3520,3615,3838,93,572,1179,
1914,2777,3768,791,2038,3413,820,2451,114,2001,
4016,2063,238,2637,1068,3723,2410,1225,168,3335,
2534,1861,1316,899,610,449,416,511,734,1085,
1564,2171,2906,3769,664,1783,3030,309,1812,3443,
1106,2993,912,3055,1230,3629,2060,619,3402,2217,
1160,231,3526,2853,2308,1891,1602,1441,1408,1503,
1726,2077,2556,3163,3898,665,1656,2775,4022,1301,
2804,339,2098,3985,1904,4047,2222,525,3052,1611,
298,3209,2152,1223,422,3845,3300,2883,2594,2433,
2400,2495,2718,3069,3548,59,794,1657,2648,3767,
918,2293,3796,1331,3090,881,2896,943,3214,1517,
4044,2603,1290,105,3144,2215,1414,741,196,3875,
3586,3425,3392,3487,3710,4061,444,1051,1786,2649,
3640,663,1910,3285,692,2323,4082,1873,3888,1935,
110,2509,940,3595,2282,1097,40,3207,2406,1733,
1188,771,482,321,288,383,606,957,1436,2043,
2778,3641,536,1655,2902,181,1684,3315,978,2865,
784,2927,1102,3501,1932,491,3274,2089,1032,103,
3398,2725,2180,1763,1474,1313,1280,1375,1598,1949,
2428,3035,3770,537,1528,2647,3894,1173,2676,211,
1970,3857,1776,3919,2094,397,2924,1483,170,3081,
2024,1095,294,3717,3172,2755,2466,2305,2272,2367,
2590,2941,3420,4027,666,1529,2520,3639,790,2165,
3668,1203,2962,753,2768,815,3086,1389,3916,2475,
1162,4073,3016,2087,1286,613,68,3747,3458,3297,
3264,3359,3582,3933,316,923,1658,2521,3512,535,
1782,3157,564,2195,3954,1745,3760,1807,4078,2381,
812,3467,2154,969,4008,3079,2278,1605,1060,643,
354,193,160,255,478,829,1308,1915,2650,3513,
408,1527,2774,53,1556,3187,850,2737,656,2799,
974,3373,1804,363,3146,1961,904,4071,3270,2597,
2052,1635,1346,1185,1152,1247,1470,1821,2300,2907,
3642,409,1400,2519,3766,1045,2548,83,1842,3729,
1648,3791,1966,269,2796,1355,42,2953,1896,967,
166,3589,3044,2627,2338,2177,2144,2239,2462,2813,
3292,3899,538,1401,2392,3511,662,2037,3540,1075,
2834,625,2640,687,2958,1261,3788,2347,1034,3945,
2888,1959,1158,485,4036,3619,3330,3169,3136,3231,
3454,3805,188,795,1530,2393,3384,407,1654,3029,
436,2067,3826,1617,3632,1679,3950,2253,684,3339,
2026,841,3880,2951,2150,1477,932,515,226,65,
32,127,350,701,1180,1787,2522,3385,280,1399,
2646,4021,1428,3059,722,2609,528,2671,846,3245,
1676,235,3018,1833,776,3943,3142,2469,1924,1507,
1218,1057,1024,1119,1342,1693,2172,2779,3514,281,
1272,2391,3638,917,2420,4051,1714,3601,1520,3663,
1838,141,2668,1227,4010,2825,1768,839,38,3461,
2916,2499,2210,2049,2016,2111,2334,2685,3164,3771,
410,1273,2264,3383,534,1909,3412,947,2706,497,
2512,559,2830,1133,3660,2219,906,3817,2760,1831,
1030,357,3908,3491,3202,3041,3008,3103,3326,3677,
60,667,1402,2265,3256,279,1526,2901,308,1939,
3698,1489,3504,1551,3822,2125,556,3211,1898,713,
3752,2823,2022,1349,804,387,98,4033,4000,4095,
222,573,1052,1659,2394,3257,152,1271,2518,3893,
1300,2931,594,2481,400,2543,718,3117,1548,107,
2890,1705,648,3815,3014,2341,1796,1379,1090,929,
896,991,1214,1565,2044,2651,3386,153,1144,2263,
3510,789,2292,3923,1586,3473,1392,3535,1710,13,
2540,1099,3882,2697,1640,711,4006,3333,2788,2371,
2082,1921,1888,1983,2206,2557,3036,3643,282,1145,
2136,3255,406,1781,3284,819,2578,369,2384,431,
2702,1005,3532,2091,778,3689,2632,1703,902,229,
3780,3363,3074,2913,2880,2975,3198,3549,4028,539,
1274,2137,3128,151,1398,2773,180,1811,3570,1361,
3376,1423,3694,1997,428,3083,1770,585,3624,2695,
1894,1221,676,259,4066,3905,3872,3967,94,445,
924,1531,2266,3129,24,1143,2390,3765,1172,2803,
466,2353,272,2415,590,2989,1420,4075,2762,1577,
520,3687,2886,2213,1668,1251,962,801,768,863,
1086,1437,1916,2523,3258,25,1016,2135,3382,661,
2164,3795,1458,3345,1264,3407,1582,3981,2412,971,
3754,2569,1512,583,3878,3205,2660,2243,1954,1793,
1760,1855,2078,2429,2908,3515,154,1017,2008,3127,
278,1653,3156,691,2450,241,2256,303,2574,877,
3404,1963,650,3561,2504,1575,774,101,3652,3235,
2946,2785,2752,2847,3070,3421,3900,411,1146,2009,
3000,23,1270,2645,52,1683,3442,1233,3248,1295,
3566,1869,300,2955,1642,457,3496,2567,1766,1093,
548,131,3938,3777,3744,3839,4062,317,796,1403,
2138,3001,3992,1015,2262,3637,1044,2675,338,2225,
144,2287,462,2861,1292,3947,2634,1449,392,3559,
2758,2085,1540,1123,834,673,640,735,958,1309,
1788,2395,3130,3993,888,2007,3254,533,2036,3667,
1330,3217,1136,3279,1454,3853,2284,843,3626,2441,
1384,455,3750,3077,2532,2115,1826,1665,1632,1727,
1950,2301,2780,3387,26,889,1880,2999,150,1525,
3028,563,2322,113,2128,175,2446,749,3276,1835,
522,3433,2376,1447,646,4069,3524,3107,2818,2657,
2624,2719,2942,3293,3772,283,1018,1881,2872,3991,
1142,2517,4020,1555,3314,1105,3120,1167,3438,1741,
172,2827,1514,329,3368,2439,1638,965,420,3,
3810,3649,3616,3711,3934,189,668,1275,2010,2873,
3864,887,2134,3509,916,2547,210,2097,16,2159,
334,2733,1164,3819,2506,1321,264,3431,2630,1957,
1412,995,706,545,512,607,830,1181,1660,2267,
3002,3865,760,1879,3126,405,1908,3539,1202,3089,
1008,3151,1326,3725,2156,715,3498,2313,1256,327,
3622,2949,2404,1987,1698,1537,1504,1599,1822,2173,
2652,3259,3994,761,1752,2871,22,1397,2900,435,
2194,4081,2000,47,2318,621,3148,1707,394,3305,
2248,1319,518,3941,3396,2979,2690,2529,2496,2591,
2814,3165,3644,155,890,1753,2744,3863,1014,2389,
3892,1427,3186,977,2992,1039,3310,1613,44,2699,
1386,201,3240,2311,1510,837,292,3971,3682,3521,
3488,3583,3806,61,540,1147,1882,2745,3736,759,
2006,3381,788,2419,82,1969,3984,2031,206,2605,
1036,3691,2378,1193,136,3303,2502,1829,1284,867,
578,417,384,479,702,1053,1532,2139,2874,3737,
632,1751,2998,277,1780,3411,1074,2961,880,3023,
1198,3597,2028,587,3370,2185,1128,199,3494,2821,
2276,1859,1570,1409,1376,1471,1694,2045,2524,3131,
3866,633,1624,2743,3990,1269,2772,307,2066,3953,
1872,4015,2190,493,3020,1579,266,3177,2120,1191,
390,3813,3268,2851,2562,2401,2368,2463,2686,3037,
3516,27,762,1625,2616,3735,886,2261,3764,1299,
3058,849,2864,911,3182,1485,4012,2571,1258,73,
3112,2183,1382,709,164,3843,3554,3393,3360,3455,
3678,4029,412,1019,1754,2617,3608,631,1878,3253,
660,2291,4050,1841,3856,1903,78,2477,908,3563,
2250,1065,8,3175,2374,1701,1156,739,450,289,
256,351,574,925,1404,2011,2746,3609,504,1623,
2870,149,1652,3283,946,2833,752,2895,1070,3469,
1900,459,3242,2057,1000,71,3366,2693,2148,1731,
1442,1281,1248,1343,1566,1917,2396,3003,3738,505,
1496,2615,3862,1141,2644,179,1938,3825,1744,3887,
2062,365,2892,1451,138,3049,1992,1063,262,3685,
3140,2723,2434,2273,2240,2335,2558,2909,3388,3995,
634,1497,2488,3607,758,2133,3636,1171,2930,721,
2736,783,3054,1357,3884,2443,1130,4041,2984,2055,
1254,581,36,3715,3426,3265,3232,3327,3550,3901,
284,891,1626,2489,3480,503,1750,3125,532,2163,
3922,1713,3728,1775,4046,2349,780,3435,2122,937,
3976,3047,2246,1573,1028,611,322,161,128,223,
446,797,1276,1883,2618,3481,376,1495,2742,21,
1524,3155,818,2705,624,2767,942,3341,1772,331,
3114,1929,872,4039,3238,2565,2020,1603,1314,1153,
1120,1215,1438,1789,2268,2875,3610,377,1368,2487,
3734,1013,2516,51,1810,3697,1616,3759,1934,237,
2764,1323,10,2921,1864,935,134,3557,3012,2595,
2306,2145,2112,2207,2430,2781,3260,3867,506,1369,
2360,3479,630,2005,3508,1043,2802,593,2608,655,
2926,1229,3756,2315,1002,3913,2856,1927,1126,453,
4004,3587,3298,3137,3104,3199,3422,3773,156,763,
1498,2361,3352,375,1622,2997,404,2035,3794,1585,
3600,1647,3918,2221,652,3307,1994,809,3848,2919,
2118,1445,900,483,194,33};





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


__global__ void interLeave(double * src, double * des , unsigned int * interLeaveTable ){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    des[tid] = src[interLeaveTable[tid]];
}

__global__ void deInterLeave(double * src, double * des , unsigned int * interLeaveTable ){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;
    des[interLeaveTable[tid]] = src[tid];
}

__global__ void gammaAlpha(double * msg ,double * parity, double * L_a, double (*gamma)[8][8], BYTE (*lastState)[8],char (*lastOut)[8] ){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;

    unsigned int s0,s2;
    for (s0=0;s0<NSTATE;s0++) {
		for (s2=0;s2<NSTATE;s2++)
			gamma[tid][s0][s2]=-INIFINITY;
		gamma[tid][s0][lastState[0][s0]]=-msg[tid]+parity[tid]*lastOut[0][s0]-log(1+exp(L_a[tid]));
		gamma[tid][s0][lastState[1][s0]]=msg[tid]+parity[tid]*lastOut[1][s0]+L_a[tid]-log(1+exp(L_a[tid]));
		//gamma[tid][s0][lastState[0][s0]]=0.5;
		//gamma[tid][s0][lastState[1][s0]]=-0.5;
    }
}

__global__ void gammaBeta(double * msg ,double * parity, double * L_a, double (*gamma)[8][8], BYTE (*nextState)[8], char (*nextOut)[8]){
    const int tid = blockIdx.x*blockDim.x + threadIdx.x;

    unsigned int s0,s2;
    for (s0=0;s0<NSTATE;s0++) {
		for (s2=0;s2<NSTATE;s2++)
			gamma[tid][s0][s2]=-INIFINITY;
		gamma[tid][s0][nextState[0][s0]]=-msg[tid]+parity[tid]*nextOut[0][s0]-log(1+exp(L_a[tid]));
		gamma[tid][s0][nextState[1][s0]]=msg[tid]+parity[tid]*nextOut[1][s0]+L_a[tid]-log(1+exp(L_a[tid]));
		//gamma[tid][s0][nextState[0][s0]]=0.5;
		//gamma[tid][s0][nextState[1][s0]]=-0.5;
    }
}

__global__ void Alpha(double (*Alpha)[8], double (*gamma)[8][8], double *maxBranch) {
	const int tid = blockIdx.x*blockDim.x + threadIdx.x;

	UINT k, s1, s2;
	double sum;

	if (tid == 0) {
		Alpha[0][0] = 0.0;
		for (s1=1;s1<NSTATE;s1++)
			Alpha[0][s1]=-INIFINITY;
	}
	else {
		for (s1=0;s1<NSTATE;s1++)
			Alpha[tid*(L_TOTAL/AlphaBetaTHREAD_NUM)][s1]=0;
	}

	//for (k=1; k<=L_TOTAL; k++) {
	for (k=tid*L_TOTAL/AlphaBetaTHREAD_NUM+1; k<(tid*L_TOTAL/AlphaBetaTHREAD_NUM+L_TOTAL/AlphaBetaTHREAD_NUM); k++) {
        for (s2=0;s2<NSTATE;s2++){
            sum = 0.0;
            for (s1=0;s1<NSTATE;s1++) {
                sum+=exp(gamma[k-1][s2][s1]+Alpha[k-1][s1]);
			}
            if (sum<MIN)
            //if (sum<=0.000000000000000000000000000001)
                Alpha[k][s2]=-INIFINITY;
            else
                Alpha[k][s2]=log(sum);
        }

		// normalization,prevent overflow
		maxBranch[k]=Alpha[k][0];
		for (s2=1;s2<NSTATE;s2++)
			if (Alpha[k][s2]>maxBranch[k])
				maxBranch[k]=Alpha[k][s2];

		for (s2=0;s2<NSTATE;s2++)
			Alpha[k][s2]=Alpha[k][s2]-maxBranch[k];
	}

}

__global__ void Beta(double (*Beta)[8], double (*gamma)[8][8], bool index, double* maxBranch) {
	const int tid = blockIdx.x*blockDim.x + threadIdx.x;

	UINT k, s1, s2;
	double sum;

	if (tid == (AlphaBetaTHREAD_NUM-1)) {
		if (index){// true -- terminated,false -- open
        Beta[L_TOTAL][0]=0.0;
        for (s2=1;s2<NSTATE;s2++)
            Beta[L_TOTAL][s2]=-INIFINITY;
		}
		else 
			for (s2=0;s2<NSTATE;s2++)
				Beta[L_TOTAL][s2]=0;
	}
	else {
		for (s2=0; s2<NSTATE; s2++)
			Beta[(tid+1)*L_TOTAL/AlphaBetaTHREAD_NUM][s2]=0;
	}

    for (k=(tid+1)*L_TOTAL/AlphaBetaTHREAD_NUM-1;k>(tid*L_TOTAL/AlphaBetaTHREAD_NUM);k--) {
   // for (k=L_TOTAL-1;k>0;k--) {

        for (s1=0;s1<NSTATE;s1++) {
            sum = 0.0;
            for (s2=0;s2<NSTATE;s2++) 
                sum += exp(gamma[k][s1][s2] + Beta[k+1][s2]);
            if (sum<MIN)
            //if (sum<=0.000000000000000000000000000001)
                Beta[k][s1] = -INIFINITY;
            else 
                Beta[k][s1] = log(sum);
        }

		// normalization,prevent overflow
		for (s2=0;s2<NSTATE;s2++)
			Beta[k][s2]=Beta[k][s2]-maxBranch[k];
	}
}

void computeAlpha(double (*AlphaHost)[8], double (*gamma)[8][8], double *maxBranch) {
    // initialize Alpha & Beta
    AlphaHost[0][0]=0;
	UINT s1,k,s2;
	double sum;
    for (s1=1;s1<NSTATE;s1++)
        AlphaHost[0][s1]=-INIFINITY;

    for (k=1;k<=L_TOTAL;k++){

        for (s2=0;s2<NSTATE;s2++){
            sum = 0;
            for (s1=0;s1<NSTATE;s1++) {
                sum+=exp(gamma[k-1][s2][s1]+AlphaHost[k-1][s1]);
			}
            if (sum<MIN)
            //if (sum<=0.000000000000000000000000000001)
                AlphaHost[k][s2]=-INIFINITY;
            else
                AlphaHost[k][s2]=log(sum);
        }

		// normalization,prevent overflow
		maxBranch[k]=AlphaHost[k][0];
		for (s2=1;s2<NSTATE;s2++)
			if (AlphaHost[k][s2]>maxBranch[k])
				maxBranch[k]=AlphaHost[k][s2];

		for (s2=0;s2<NSTATE;s2++)
			AlphaHost[k][s2]=AlphaHost[k][s2]-maxBranch[k];
    }
}

void computeBeta(double (*BetaHost)[8], double (*gamma)[8][8], bool index, double * maxBranch){
    // initialize Beta
	UINT s1,k,s2;
	double sum;
    if (index){// true -- terminated,false -- open
        BetaHost[L_TOTAL][0]=0;
        for (s2=1;s2<NSTATE;s2++)
            BetaHost[L_TOTAL][s2]=-INIFINITY;
    }
    else 
        for (s2=0;s2<NSTATE;s2++)
            BetaHost[L_TOTAL][s2]=0;

    for (k=L_TOTAL-1;k>0;k--) {

        for (s1=0;s1<NSTATE;s1++) {
            sum = 0.0;
            for (s2=0;s2<NSTATE;s2++) 
                sum += exp(gamma[k][s1][s2] + BetaHost[k+1][s2]);
            if (sum<MIN)
            //if (sum<=0.000000000000000000000000000001)
                BetaHost[k][s1] = -INIFINITY;
            else 
                BetaHost[k][s1] = log(sum);
        }

		// normalization,prevent overflow
		for (s2=0;s2<NSTATE;s2++)
			BetaHost[k][s2]=BetaHost[k][s2]-maxBranch[k];
    }
}

__global__ void normalizationAlphaAndBeta(double (*Alpha)[8], double (*Beta)[8]) {
    unsigned int tid = threadIdx.x+1; 
    double max_branch;
    max_branch = Alpha[tid][0];
	UINT s2;
    for (s2=1;s2<NSTATE;s2++)
        if (Alpha[tid][s2]>max_branch)
            max_branch = Alpha[tid][s2];

    for (s2=0;s2<NSTATE;s2++) {
        Alpha[tid][s2] = Alpha[tid][s2] - max_branch;

        if (tid != L_TOTAL) 
            Beta[tid][s2] = Beta[tid][s2] - max_branch;
    }

}

__global__ void LLRS(double * msg, double * parity, double * L_a, double (*Alpha)[8], double (*Beta)[8], double * L_all,BYTE (*lastState)[8], char (*lastOut)[8]) {
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    UINT s2;
	double sum0 = 0.0, sum1 = 0.0;
    for (s2=0;s2<NSTATE;s2++) {
        //gamma[LastState[0][s2]]=-msg[tid]+parity[tid]*LastOut[0][s2]-log(1+exp(L_a[tid]));
        //gamma[LastState[1][s2]]=msg[tid]+parity[tid]*LastOut[1][s2]+L_a[tid]-log(1+exp(L_a[tid]));
        //sum0+=exp(gamma[LastState[0][s2]]+Alpha[tid][LastState[0][s2]]+Beta[tid+1][s2]);
        //sum1+=exp(gamma[LastState[1][s2]]+Alpha[tid][LastState[1][s2]]+Beta[tid+1][s2]);
        double gamma0=-msg[tid]+parity[tid]*lastOut[0][s2]-log(1+exp(L_a[tid]));
        double gamma1=msg[tid]+parity[tid]*lastOut[1][s2]+L_a[tid]-log(1+exp(L_a[tid]));
        sum0+=exp(gamma0+Alpha[tid][lastState[0][s2]]+Beta[tid+1][s2]);
        sum1+=exp(gamma1+Alpha[tid][lastState[1][s2]]+Beta[tid+1][s2]);
    }
    //L_all[tid]=log(sum1)-log(sum0);
    L_all[tid]=log(sum1)-log(sum0);
}

__global__ void extrinsicInformation(double * L_all, double * msg, double * L_a, double * L_e) {
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    L_e[tid] = L_all[tid] - 2*msg[tid] - L_a[tid];
}

__global__ void demultiplex(double * stream, double * msg, double * parity0, double * parity1) {
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
        msg[tid]=stream[3*tid];
        parity0[tid]=stream[3*tid+1];
        parity1[tid]=stream[3*tid+2];
}

__global__ void initializeExtrinsicInformation(double * L_e) {
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    L_e[tid] = 0;
    
}

__global__ void exestimateInformationBits(double * L_all, BYTE * msghat, UINT * m_Inter_table) {
    unsigned int tid = blockIdx.x*blockDim.x + threadIdx.x;
    if(L_all[tid]>0)
        msghat[m_Inter_table[tid]]=1;
    else
        msghat[m_Inter_table[tid]]=0;
}

void countErrors(BYTE *m, BYTE * mhat, UINT * bitsError, UINT * frameError, UINT iter) {

	bool f_err = false;
	for (int i=0; i<(L_TOTAL-M);i++) {
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
    initRandom(seed);

	BYTE * m;
	BYTE * x;
	double * y;
	BYTE * mhat;

	int frame;
	UINT bits_all,bits_err[MAXITER],frame_err[MAXITER];
	double Ber,Fer;
	double Eb_No_dB,No;
	//bool f_err;
	//FILE * fp;
	int i;

	m = new BYTE[L_TOTAL];
	x = new BYTE[L_ALL];
	y = new double[L_ALL];
	mhat = new BYTE[L_TOTAL];

	//init_Block_interleave_table();	// block interleave


    
    findCudaDevice(argc, (const char **)argv);

	BYTE (*LastStateDevice)[8];
	BYTE (*NextStateDevice)[8];
	char (*LastOutDevice)[8];
	char (*NextOutDevice)[8];
	double * yDevice;
	double * msgDevice;
	double * imsgDevice;
	BYTE * mhatDevice;
	double * parity0Device;
	double * parity1Device;
	UINT * tableDevice;
	double * L_eDevice;
	double * L_aDevice;
	double * L_allDevice;
	double (*gammaAlphaDevice)[8][8];
	double (*gammaBetaDevice)[8][8];
	double (*AlphaDevice)[8];
	double (*BetaDevice)[8];
	double *maxBranchDevice;

    cudaMalloc((void **)&LastStateDevice, 2*8*sizeof(BYTE));
    cudaMalloc((void **)&NextStateDevice, 2*8*sizeof(BYTE));
    cudaMalloc((void **)&LastOutDevice, 2*8*sizeof(char));
    cudaMalloc((void **)&NextOutDevice, 2*8*sizeof(char));

    cudaMalloc((void **)&yDevice, L_ALL*sizeof(double));
    cudaMalloc((void **)&msgDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&imsgDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&mhatDevice, L_TOTAL*sizeof(BYTE));
    cudaMalloc((void **)&parity0Device, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&parity1Device, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&tableDevice, L_TOTAL*sizeof(unsigned int));
    cudaMalloc((void **)&L_eDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&L_aDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&L_allDevice, L_TOTAL*sizeof(double));
    cudaMalloc((void **)&gammaAlphaDevice, L_TOTAL*sizeof(double)*8*8);
    cudaMalloc((void **)&gammaBetaDevice, L_TOTAL*sizeof(double)*8*8);
    cudaMalloc((void **)&AlphaDevice, (L_TOTAL+1)*sizeof(double)*8);
    cudaMalloc((void **)&BetaDevice, (L_TOTAL+1)*sizeof(double)*8);
    cudaMalloc((void **)&maxBranchDevice, (L_TOTAL+1)*sizeof(double));

	//For Debug
	//double L_aHost[L_TOTAL];
	//double L_aHost1[L_TOTAL];
	//double L_allHost[L_TOTAL];

    double gammaAlphaHost[L_TOTAL][8][8];
    double gammaBetaHost[L_TOTAL][8][8];
    double AlphaHost[L_TOTAL+1][8];
    double BetaHost[L_TOTAL+1][8];

	double max_branch[L_TOTAL+1];

    cudaMemcpy(LastStateDevice,LastState,sizeof(BYTE)*2*8, cudaMemcpyHostToDevice);
    cudaMemcpy(NextStateDevice,NextState,sizeof(BYTE)*2*8, cudaMemcpyHostToDevice);
    cudaMemcpy(LastOutDevice,LastOut,sizeof(char)*2*8, cudaMemcpyHostToDevice);
    cudaMemcpy(NextOutDevice,NextOut,sizeof(char)*2*8, cudaMemcpyHostToDevice);

    cudaMemcpy(tableDevice,m_Inter_table,sizeof(unsigned int)*L_TOTAL, cudaMemcpyHostToDevice);

	for (Eb_No_dB= 0.0;Eb_No_dB<3.0;Eb_No_dB+=0.5){

	//Eb_No_dB = 0.0;
		No = 1/pow(10.0,Eb_No_dB/10.0);
		bits_all = 0;
		for (i =0; i<MAXITER;i++) {
			bits_err[i]=0;
			frame_err[i]=0;
		}

		for (frame = 0; frame<FRAME_NUM; frame++, bits_all += (L_TOTAL-M)) {

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

			cudaMemcpy(yDevice,y,sizeof(double)*L_ALL, cudaMemcpyHostToDevice);

			demultiplex<<<BLOCK_NUM,THREAD_NUM>>>(yDevice, msgDevice, parity0Device, parity1Device); 
			interLeave<<<BLOCK_NUM,THREAD_NUM>>>(msgDevice, imsgDevice, tableDevice);
			initializeExtrinsicInformation<<<BLOCK_NUM,THREAD_NUM>>>(L_eDevice);

			for (int iter = 0; iter<MAXITER; iter++) {
				
				deInterLeave<<<BLOCK_NUM,THREAD_NUM>>>(L_eDevice, L_aDevice, tableDevice);

				gammaAlpha<<<BLOCK_NUM,THREAD_NUM>>>(msgDevice , parity0Device,  L_aDevice,  gammaAlphaDevice,LastStateDevice, LastOutDevice);
				gammaBeta<<<BLOCK_NUM,THREAD_NUM>>>(msgDevice , parity0Device,  L_aDevice,  gammaBetaDevice, NextStateDevice, NextOutDevice);
				Alpha<<<1,AlphaBetaTHREAD_NUM>>>(AlphaDevice, gammaAlphaDevice, maxBranchDevice);
				Beta<<<1,AlphaBetaTHREAD_NUM>>>(BetaDevice, gammaBetaDevice,true, maxBranchDevice);
				//cudaMemcpy(AlphaHost, AlphaDevice, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyDeviceToHost);
				//cudaMemcpy(gammaAlphaHost, gammaAlphaDevice, sizeof(double)*L_TOTAL*8*8, cudaMemcpyDeviceToHost);
				//cudaMemcpy(gammaBetaHost, gammaBetaDevice, sizeof(double)*L_TOTAL*8*8, cudaMemcpyDeviceToHost);

				//computeAlpha(AlphaHost, gammaAlphaHost, max_branch);
				//computeBeta(BetaHost, gammaBetaHost, true,max_branch);
				//cudaMemcpy(AlphaDevice, AlphaHost, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyHostToDevice);
				//cudaMemcpy(BetaDevice, BetaHost, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyHostToDevice);
				//normalizationAlphaAndBeta<<<BLOCK_NUM,THREAD_NUM>>>(AlphaDevice, BetaDevice);
				//cudaMemcpy(AlphaHost, AlphaDevice, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyDeviceToHost);

				LLRS<<<BLOCK_NUM,THREAD_NUM>>>(msgDevice, parity0Device, L_aDevice, AlphaDevice, BetaDevice, L_allDevice,LastStateDevice, LastOutDevice);

				extrinsicInformation<<<BLOCK_NUM,THREAD_NUM>>>(L_allDevice, msgDevice, L_aDevice, L_eDevice);
				//if (iter >= 3) {
				///debug
				//exestimateInformationBits<<<BLOCK_NUM,THREAD_NUM>>>(L_allDevice, mhatDevice, tableDevice); 

				//cudaMemcpy(mhat, mhatDevice, sizeof(BYTE)*L_TOTAL, cudaMemcpyDeviceToHost);
				//countErrors(m, mhat, bits_err, frame_err, iter);
				//debug
				//cudaMemcpy(L_aHost1, L_aDevice, sizeof(double)*L_TOTAL, cudaMemcpyDeviceToHost);
				//}

				interLeave<<<BLOCK_NUM,THREAD_NUM>>>(L_eDevice, L_aDevice, tableDevice);

				gammaAlpha<<<BLOCK_NUM,THREAD_NUM>>>(imsgDevice , parity1Device,  L_aDevice,  gammaAlphaDevice, LastStateDevice, LastOutDevice);
				gammaBeta<<<BLOCK_NUM,THREAD_NUM>>>(imsgDevice , parity1Device,  L_aDevice,  gammaBetaDevice, NextStateDevice, NextOutDevice);
				Alpha<<<1,AlphaBetaTHREAD_NUM>>>(AlphaDevice, gammaAlphaDevice, maxBranchDevice);
				Beta<<<1,AlphaBetaTHREAD_NUM>>>(BetaDevice, gammaBetaDevice,false, maxBranchDevice);
				//cudaMemcpy(AlphaHost, AlphaDevice, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyDeviceToHost);
				//cudaMemcpy(gammaAlphaHost, gammaAlphaDevice, sizeof(double)*L_TOTAL*8*8, cudaMemcpyDeviceToHost);
				//cudaMemcpy(gammaBetaHost, gammaBetaDevice, sizeof(double)*L_TOTAL*8*8, cudaMemcpyDeviceToHost);

				//computeAlpha(AlphaHost, gammaAlphaHost, max_branch);
				//Beta(BetaHost, gammaBetaHost, false, max_branch);
				//cudaMemcpy(AlphaDevice, AlphaHost, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyHostToDevice);
				//cudaMemcpy(BetaDevice, BetaHost, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyHostToDevice);
				//normalizationAlphaAndBeta<<<BLOCK_NUM,THREAD_NUM>>>(AlphaDevice, BetaDevice);
				//cudaMemcpy(AlphaHost, AlphaDevice, sizeof(double)*(L_TOTAL+1)*8, cudaMemcpyDeviceToHost);

				LLRS<<<BLOCK_NUM,THREAD_NUM>>>(imsgDevice, parity1Device, L_aDevice, AlphaDevice, BetaDevice, L_allDevice, LastStateDevice,LastOutDevice);

				extrinsicInformation<<<BLOCK_NUM,THREAD_NUM>>>(L_allDevice, imsgDevice, L_aDevice, L_eDevice);

				exestimateInformationBits<<<BLOCK_NUM,THREAD_NUM>>>(L_allDevice, mhatDevice, tableDevice); 

				cudaMemcpy(mhat, mhatDevice, sizeof(BYTE)*L_TOTAL, cudaMemcpyDeviceToHost);
				countErrors(m, mhat, bits_err, frame_err, iter);

				//debug
				//cudaMemcpy(L_aHost, L_aDevice, sizeof(double)*L_TOTAL, cudaMemcpyDeviceToHost);
			}
			// estimate information bits
			//exestimateInformationBits<<<BLOCK_NUM,THREAD_NUM>>>(L_allDevice, mhatDevice, tableDevice); 

			//cudaMemcpy(mhat, mhatDevice, sizeof(BYTE)*L_TOTAL, cudaMemcpyDeviceToHost);
			// count errors
			//UINT bits_err = 0;
			//for (i=0;i<L_TOTAL-M;i++) {
			//	if (mhat[i]!=m[i]) {
			//		bits_err++;
			//	}
			//}
			//cout<<"bits_err: "<<bits_err<<endl;
		}

		printf("-------------------------\n");
		printf("Eb/No=%fdB:\n",Eb_No_dB);
		printf("-------------------------\n");
		//fprintf(fp,"-------------------------\n");
		//fprintf(fp,"Eb/No=%fdB:\n",Eb_No_dB);
		//fprintf(fp,"-------------------------\n");

		for (i=0;i<MAXITER;i++) {
			Ber=(double)bits_err[i]/(double)bits_all;
			Fer=(double)frame_err[i]/(double)FRAME_NUM;
			printf("Iteration:%d\n",i+1);
			printf("---Ber=%f\n---Fer=%f\n",Ber,Fer);
			//fprintf(fp,"Iteration:%d\n",i);
			//fprintf(fp,"---Ber=%f\n---Fer=%f\n",Ber,Fer);
		}
	}

	delete m;
	delete x;
	delete y;
	delete mhat;

	cudaFree(LastStateDevice);
	cudaFree(NextStateDevice);
	cudaFree(LastOutDevice);
	cudaFree(NextOutDevice);
	cudaFree(yDevice);
	cudaFree(msgDevice);
	cudaFree(imsgDevice);
	cudaFree(mhatDevice);
	cudaFree(parity0Device);
	cudaFree(parity1Device);
	cudaFree(tableDevice);
	cudaFree(L_eDevice);
	cudaFree(L_aDevice);
	cudaFree(L_allDevice);
	cudaFree(gammaAlphaDevice);
	cudaFree(gammaBetaDevice);
	cudaFree(AlphaDevice);
	cudaFree(BetaDevice);

	return 0;
}
