#pragma once
#include <stdint.h>


// Functions optimized for worst-case error


//Searching constants between 0x1fba7da0 and 0x1fc00000...
//  ...exhaustive between 1fbb66c4 and 1fbb68c5...
//  ...best constant is 1fbb67a5 with error score 0.000601071
/*
        approximate x^0.5 with 1 newtonian steps, x=[1,4]
                RMS error: 0.000228914
                mean error: 0.000179662
                worst error: 0.000601071
*/
float rb_2_root(const float y)
{
        union {float x; int32_t i;}; x = y; // interpret float as integer
        i = 0x1fbb67a5 + (i >> 1); // log-approximation hack
        x = 0.5f * x + 0.5f * y / x; // newtonian step #1
        return x;
}

//Searching constants between 0x5f2f7900 and 0x5f400000...
//  ...exhaustive between 5f375855 and 5f375a56...
//  ...best constant is 5f375a55 with error score 0.00175157
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

//Searching constants between 0x2a4dfcc0 and 0x2a555540...
//  ...exhaustive between 2a511f6c and 2a51216e...
//  ...best constant is 2a51206c with error score 0.000993095
/*
        approximate x^0.333333 with 1 newtonian steps, x=[1,8]
                RMS error: 0.000410751
                mean error: 0.000325521
                worst error: 0.000993095
*/
float rb_3_root(const float y)
{
        union {float x; int32_t i;}; x = y; // interpret float as integer
        i = 0x2a51206c + (i / 3); // log-approximation hack
        x = 0.666667f * x + 0.333333f * y / (x*x); // newtonian step #1
        return x;
}

//Searching constants between 0x549bfa00 and 0x54aaab00...
//  ...exhaustive between 54a21cfa and 54a21efb...
//  ...best constant is 54a21e32 with error score 0.00233629
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

//Searching constants between 0x2f97bc80 and 0x2fa00000...
//  ...exhaustive between 2f9bdd40 and 2f9bdf41...
//  ...best constant is 2f9bdd40 with error score 0.00201698
/*
        approximate x^0.25 with 1 newtonian steps, x=[1,16]
                RMS error: 0.000679
                mean error: 0.000488281
                worst error: 0.00201698
*/
float rb_4_root(const float y)
{
        union {float x; int32_t i;}; x = y; // interpret float as integer
        i = 0x2f9bdd40 + (i >> 2); // log-approximation hack
        x = 0.75f * x + 0.25f * y / (x*x*x); // newtonian step #1
        return x;
}

//Searching constants between 0x4f523a00 and 0x4f600000...
//  ...exhaustive between 4f583f9e and 4f5841a0...
//  ...best constant is 4f5841a0 with error score 0.00243795
/*
        approximate x^-0.25 with 1 newtonian steps, x=[1,16]
                RMS error: 0.00129543
                mean error: -0.000951638
                worst error: -0.00243795
*/
float rb_inv_4_root(const float y)
{
        union {float x; int32_t i;}; x = y; // interpret float as integer
        i = 0x4f5841a0 - (i >> 2); // log-approximation hack
        x *= 1.25f - 0.25f * y * (x*x*x*x); // newtonian step #1
        return x;
}