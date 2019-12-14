#pragma once


// Functions optimized for worst-case error

/*
	approximate x^0.5 with 1 newtonian steps, x=[1,4]
		RMS error: 0.000228909
		mean error: 0.00017966
		worst error: 0.000601098
*/
float rb_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x1fbb67a9 + (i >> 1); // log-approximation hack
	x *= 0.5f + 0.5f * y / (x*x); // newtonian step #1
	return x;
}

/*
	approximate x^-0.5 with 1 newtonian steps, x=[1,4]
		RMS error: 0.00108793
		mean error: -0.000980849
		worst error: -0.00175157
*/
float rb_inv_2_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x5f375a55 - (i >> 1); // log-approximation hack
	x *= 1.5f - 0.5f * y * (x*x); // newtonian step #1
	return x;
}

/*
	approximate x^0.333333 with 1 newtonian steps, x=[1,8]
		RMS error: 0.000410755
		mean error: 0.000325521
		worst error: 0.000993097
*/
float rb_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2a512072 + (i / 3); // log-approximation hack
	x *= 0.666667f + 0.333333f * y / (x*x*x); // newtonian step #1
	return x;
}

/*
	approximate x^-0.333333 with 1 newtonian steps, x=[1,8]
		RMS error: 0.00114492
		mean error: -0.000763734
		worst error: -0.00233629
*/
float rb_inv_3_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x54a21e32 - (i / 3); // log-approximation hack
	x *= 1.33333f - 0.333333f * y * (x*x*x); // newtonian step #1
	return x;
}

/*
	approximate x^0.25 with 1 newtonian steps, x=[1,16]
		RMS error: 0.000679
		mean error: 0.000488281
		worst error: 0.0020169
*/
float rb_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x2f9bdd40 + (i >> 2); // log-approximation hack
	x *= 0.75f + 0.25f * y / (x*x*x*x); // newtonian step #1
	return x;
}

/*
	approximate x^-0.25 with 1 newtonian steps, x=[1,16]
		RMS error: 0.00129543
		mean error: -0.000951639
		worst error: -0.00243795
*/
float rb_inv_4_root(const float y)
{
	union {float x; int32_t i;}; x = y; // interpret float as integer
	i = 0x4f5841a0 - (i >> 2); // log-approximation hack
	x *= 1.25f - 0.25f * y * (x*x*x*x); // newtonian step #1
	return x;
}
