#pragma once
#include <stdint.h>


// Functions optimized for worst-case error


//Searching k in [0x1fba7da0,0x1fc00000], m in [0.25,1]...
//  ...best design k=1fbed46c, m=0.510934 with error score 0.000239134
/*
	Approximate x^(1/2) with 1 newtonian steps
	Error:
		RMS:  0.000148086
		mean: -4.55815e-05
		min:  -0.000239272 @ 1.62038
		max:  0.000239223 @ 1.01828
*/
float rb_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x1fbed46c + (i >> 1); // log-approximation hack
	x = 0.489066f * x + 0.510934f * y / x; // newtonian step #1
	return x;
}

//Searching k in [0x5f2f7900,0x5f400000], m in [-0.25,-1]...
//  ...best design k=5f329cb4, m=-0.535152 with error score 0.000775337
/*
	Approximate x^(1/-2) with 1 newtonian steps
	Error:
		RMS:  0.000492882
		mean: -2.13658e-05
		min:  -0.000775413 @ 3.58168
		max:  0.000775621 @ 3.22598
*/
float rb_inv_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x5f329cb4 - (i >> 1); // log-approximation hack
	x *= 1.53515f - 0.535152f * y * (x*x); // newtonian step #1
	return x;
}

//Searching k in [0x2a4dfcc0,0x2a555540], m in [0.166667,0.666667]...
//  ...best design k=2a543a8a, m=0.347254 with error score 0.000430107
/*
	Approximate x^(1/3) with 1 newtonian steps
	Error:
		RMS:  0.000237942
		mean: -9.79593e-05
		min:  -0.000430293 @ 1.63116
		max:  0.000430241 @ 1.02589
*/
float rb_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2a543a8a + (i / 3); // log-approximation hack
	x = 0.652746f * x + 0.347254f * y / (x*x); // newtonian step #1
	return x;
}

//Searching k in [0x549bfa00,0x54aaab00], m in [-0.166667,-0.666667]...
//  ...best design k=549da829, m=-0.364694 with error score 0.00102782
/*
	Approximate x^(1/-3) with 1 newtonian steps
	Error:
		RMS:  0.000742797
		mean: 0.000277367
		min:  -0.00102797 @ 6.78041
		max:  0.00102644 @ 5.40779
*/
float rb_inv_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x549da829 - (i / 3); // log-approximation hack
	x *= 1.36469f - 0.364694f * y * (x*x*x); // newtonian step #1
	return x;
}

//Searching k in [0x2f97bc80,0x2fa00000], m in [0.125,0.5]...
//  ...best design k=2f9ed7c2, m=0.266598 with error score 0.000713944
/*
	Approximate x^(1/4) with 1 newtonian steps
	Error:
		RMS:  0.000444153
		mean: -0.000195166
		min:  -0.000714091 @ 1.71197
		max:  0.000714014 @ 1.03616
*/
float rb_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2f9ed7c2 + (i >> 2); // log-approximation hack
	x = 0.733402f * x + 0.266598f * y / (x*x*x); // newtonian step #1
	return x;
}

//Searching k in [0x4f523a00,0x4f600000], m in [-0.125,-0.5]...
//  ...best design k=4f541f5e, m=-0.277463 with error score 0.0011096
/*
	Approximate x^(1/-4) with 1 newtonian steps
	Error:
		RMS:  0.000731492
		mean: 0.000177324
		min:  -0.0011097 @ 13.0308
		max:  0.00110985 @ 1.32538
*/
float rb_inv_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x4f541f5e - (i >> 2); // log-approximation hack
	x *= 1.27746f - 0.277463f * y * (x*x*x*x); // newtonian step #1
	return x;
}

