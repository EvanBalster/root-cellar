#pragma once
#include <stdint.h>


// Functions optimized for worst-case error


//Searching k in [0x1fba7da0,0x1fc00000], m in [0.5,0.75] ...........
//  ...best design k=1fbed49a, m=0.510929 with error score 0.000239058
/*
	Approximate x^(1/2) with 1 newtonian steps
	Error:
		RMS:  0.000148007
		mean: -4.54788e-05
		min:  -0.000239056 @ 1.62077
		max:  0.000239058 @ 1.99997
*/
float rb_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x1fbed49a + (i >> 1); // log-approximation hack
	x = 0.489071f * x + 0.510929f * y / x; // newtonian step #1
	return x;
}

//Searching k in [0x5f2f7900,0x5f400000], m in [-0.5,-0.75] ...........
//  ...best design k=5f32a121, m=-0.535102 with error score 0.000773445
/*
	Approximate x^(1/-2) with 1 newtonian steps
	Error:
		RMS:  0.000494072
		mean: -2.5526e-05
		min:  -0.000773442 @ 3.58223
		max:  0.000773445 @ 3.22659
*/
float rb_inv_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x5f32a121 - (i >> 1); // log-approximation hack
	x *= 1.5351f - 0.535102f * y * (x*x); // newtonian step #1
	return x;
}

//Searching k in [0x2a4dfcc0,0x2a555540], m in [0.333333,0.5] ...........
//  ...best design k=2a543aa3, m=0.347252 with error score 0.000430098
/*
	Approximate x^(1/3) with 1 newtonian steps
	Error:
		RMS:  0.000237859
		mean: -9.77778e-05
		min:  -0.000430098 @ 1.63107
		max:  0.000430033 @ 1.99999
*/
float rb_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2a543aa3 + (i / 3); // log-approximation hack
	x = 0.652748f * x + 0.347252f * y / (x*x); // newtonian step #1
	return x;
}

//Searching k in [0x549bfa00,0x54aaab00], m in [-0.333333,-0.5] ...........
//  ...best design k=549da7bf, m=-0.364707 with error score 0.00102717
/*
	Approximate x^(1/-3) with 1 newtonian steps
	Error:
		RMS:  0.000742809
		mean: 0.000277928
		min:  -0.00102717 @ 6.78018
		max:  0.00102716 @ 5.40931
*/
float rb_inv_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x549da7bf - (i / 3); // log-approximation hack
	x *= 1.36471f - 0.364707f * y * (x*x*x); // newtonian step #1
	return x;
}

//Searching k in [0x2f97bc80,0x2fa00000], m in [0.25,0.375] ...........
//  ...best design k=2f9ed7c0, m=0.266598 with error score 0.000714053
/*
	Approximate x^(1/4) with 1 newtonian steps
	Error:
		RMS:  0.000444122
		mean: -0.000195135
		min:  -0.000714043 @ 1.71143
		max:  0.000714053 @ 3.99998
*/
float rb_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2f9ed7c0 + (i >> 2); // log-approximation hack
	x = 0.733402f * x + 0.266598f * y / (x*x*x); // newtonian step #1
	return x;
}

//Searching k in [0x4f523a00,0x4f600000], m in [-0.25,-0.375] ...........
//  ...best design k=4f542107, m=-0.277446 with error score 0.00110848
/*
	Approximate x^(1/-4) with 1 newtonian steps
	Error:
		RMS:  0.000733642
		mean: 0.000175304
		min:  -0.00110848 @ 13.0323
		max:  0.00110847 @ 1.3249
*/
float rb_inv_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x4f542107 - (i >> 2); // log-approximation hack
	x *= 1.27745f - 0.277446f * y * (x*x*x*x); // newtonian step #1
	return x;
}


//Searching k in [0x1fba7da0,0x1fc00000], m in [0.5,0.5] .........
//  ...best design k=1fbb4f2e, m=0.5 with error score 0.0347475
/*
	Approximate x^(1/2) with 0 newtonian steps
	Error:
		RMS:  0.0190506
		mean: -0.005133
		min:  -0.0347474 @ 1.07329
		max:  0.0347475 @ 2
*/
float rb0_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x1fbb4f2e + (i >> 1); // log-approximation hack
	return x;
}

//Searching k in [0x5f2f7900,0x5f400000], m in [-0.5,-0.5] ..........
//  ...best design k=5f37642f, m=-0.5 with error score 0.0342129
/*
	Approximate x^(1/-2) with 0 newtonian steps
	Error:
		RMS:  0.0244769
		mean: 0.01338
		min:  -0.0342129 @ 3.73098
		max:  0.0342129 @ 2.57689
*/
float rb0_inv_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x5f37642f - (i >> 1); // log-approximation hack
	return x;
}

//Searching k in [0x2a4dfcc0,0x2a555540], m in [0.333333,0.333333] .........
//  ...best design k=2a510680, m=0.333333 with error score 0.0315547
/*
	Approximate x^(1/3) with 0 newtonian steps
	Error:
		RMS:  0.0180422
		mean: 0.00321401
		min:  -0.0315546 @ 1.10097
		max:  0.0315547 @ 2
*/
float rb0_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2a510680 + (i / 3); // log-approximation hack
	return x;
}

//Searching k in [0x549bfa00,0x54aaab00], m in [-0.333333,-0.333333] ..........
//  ...best design k=54a232a3, m=-0.333333 with error score 0.0342405
/*
	Approximate x^(1/-3) with 0 newtonian steps
	Error:
		RMS:  0.0195931
		mean: 0.0068072
		min:  -0.0342405 @ 7.20605
		max:  0.0342405 @ 2.90059
*/
float rb0_inv_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x54a232a3 - (i / 3); // log-approximation hack
	return x;
}

//Searching k in [0x2f97bc80,0x2fa00000], m in [0.25,0.25] .........
//  ...best design k=2f9b374e, m=0.25 with error score 0.0342323
/*
	Approximate x^(1/4) with 0 newtonian steps
	Error:
		RMS:  0.015625
		mean: 0.00425431
		min:  -0.034232 @ 1.1495
		max:  0.0342323 @ 4
*/
float rb0_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2f9b374e + (i >> 2); // log-approximation hack
	return x;
}

//Searching k in [0x4f523a00,0x4f600000], m in [-0.25,-0.25] ..........
//  ...best design k=4f58605b, m=-0.25 with error score 0.0312108
/*
	Approximate x^(1/-4) with 0 newtonian steps
	Error:
		RMS:  0.020528
		mean: 0.00879536
		min:  -0.0312108 @ 14.0941
		max:  0.0312107 @ 5.40988
*/
float rb0_inv_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x4f58605b - (i >> 2); // log-approximation hack
	return x;
}



//Searching k in [0x1fba7da0,0x1fc00000], m in [0.5,0.75] ...........
//  ...best design k=1fbb75ad, m=0.500122 with error score 1.68567e-07
/*
	Approximate x^(1/2) with 2 newtonian steps
	Error:
		RMS:  4.7967e-08
		mean: -1.41369e-08
		min:  -1.64985e-07 @ 2.08829
		max:  1.68567e-07 @ 2.00049
*/
float rb2_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x1fbb75ad + (i >> 1); // log-approximation hack
	x = 0.499878f * x + 0.500122f * y / x; // newtonian step #1
	x = 0.499878f * x + 0.500122f * y / x; // newtonian step #2
	return x;
}

//Searching k in [0x5f2f7900,0x5f400000], m in [-0.5,-0.75] ...........
//  ...best design k=5f3634f9, m=-0.501326 with error score 1.40452e-06
/*
	Approximate x^(1/-2) with 2 newtonian steps
	Error:
		RMS:  8.95917e-07
		mean: 6.21959e-07
		min:  -1.38048e-06 @ 3.7251
		max:  1.40452e-06 @ 2.83293
*/
float rb2_inv_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x5f3634f9 - (i >> 1); // log-approximation hack
	x *= 1.50133f - 0.501326f * y * (x*x); // newtonian step #1
	x *= 1.50133f - 0.501326f * y * (x*x); // newtonian step #2
	return x;
}

//Searching k in [0x2a4dfcc0,0x2a555540], m in [0.333333,0.5] ...........
//  ...best design k=2a4fcd03, m=0.333818 with error score 6.45394e-07
/*
	Approximate x^(1/3) with 2 newtonian steps
	Error:
		RMS:  2.90881e-07
		mean: -2.40255e-07
		min:  -6.45394e-07 @ 1.36116
		max:  5.73198e-07 @ 1.12441
*/
float rb2_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2a4fcd03 + (i / 3); // log-approximation hack
	x = 0.666182f * x + 0.333818f * y / (x*x); // newtonian step #1
	x = 0.666182f * x + 0.333818f * y / (x*x); // newtonian step #2
	return x;
}

//Searching k in [0x549bfa00,0x54aaab00], m in [-0.333333,-0.5] ...........
//  ...best design k=54a1b99d, m=-0.334677 with error score 2.18458e-06
/*
	Approximate x^(1/-3) with 2 newtonian steps
	Error:
		RMS:  1.04454e-06
		mean: 6.17437e-07
		min:  -2.18455e-06 @ 7.17768
		max:  2.18458e-06 @ 3.56147
*/
float rb2_inv_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x54a1b99d - (i / 3); // log-approximation hack
	x *= 1.33468f - 0.334677f * y * (x*x*x); // newtonian step #1
	x *= 1.33468f - 0.334677f * y * (x*x*x); // newtonian step #2
	return x;
}

//Searching k in [0x2f97bc80,0x2fa00000], m in [0.25,0.375] ...........
//  ...best design k=2f9b8068, m=0.250534 with error score 9.49041e-07
/*
	Approximate x^(1/4) with 2 newtonian steps
	Error:
		RMS:  5.28477e-07
		mean: -4.76837e-07
		min:  -9.49041e-07 @ 1.01967
		max:  9.27656e-07 @ 3.99267
*/
float rb2_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2f9b8068 + (i >> 2); // log-approximation hack
	x = 0.749466f * x + 0.250534f * y / (x*x*x); // newtonian step #1
	x = 0.749466f * x + 0.250534f * y / (x*x*x); // newtonian step #2
	return x;
}

//Searching k in [0x4f523a00,0x4f600000], m in [-0.25,-0.375] ...........
//  ...best design k=4f58020d, m=-0.251282 with error score 2.76944e-06
/*
	Approximate x^(1/-4) with 2 newtonian steps
	Error:
		RMS:  1.3487e-06
		mean: 5.97779e-07
		min:  -2.65373e-06 @ 14.041
		max:  2.76944e-06 @ 3.93
*/
float rb2_inv_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x4f58020d - (i >> 2); // log-approximation hack
	x *= 1.25128f - 0.251282f * y * (x*x*x*x); // newtonian step #1
	x *= 1.25128f - 0.251282f * y * (x*x*x*x); // newtonian step #2
	return x;
}
