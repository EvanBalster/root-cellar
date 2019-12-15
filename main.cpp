//
//  main.cpp
//  rootbeer
//
//  Created by Evan Balster on 19-12-13.
//  Copyright © 2019 Interactopia LLC. All rights reserved.
//

#if _MSC_VER
#include <immintrin.h>
#endif

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>

#include "root_cellar.h"
#include "root_cellar_generated.h"

using namespace rootbeer;

static float TEST_VALUES[8192];

template<typename T_Func>
void Print_Test_Root_Approx(const char *name, const T_Func &func, int root_index)
{
	using float_t = decltype(func(1.0));
	
	float_t
		exponent = float_t(1)/root_index,
		range_min = float_t(1),
		range_max = float_t(1 << std::abs(root_index));

	std::cout << "\tapproximate x^" << exponent << " with " << name
		<< ", x=[" << range_min << "," << range_max << "]" << std::endl;
	auto test = Test_Pow_Approx(func, exponent, range_min, range_max);
	std::cout
		<< "\t\tRMS error: " << std::sqrt(test.mean_sq_error) << std::endl
		<< "\t\tmean error: " << test.mean_error << std::endl
		<< "\t\tworst error: " << test.worst_error << std::endl;
}

template<typename T_Func>
void Print_Func_Profile(const char *name, const T_Func &func)
{
	float total = 0.f;
	auto start = std::chrono::high_resolution_clock::now();
	for (const auto v : TEST_VALUES)
	{
		total += func(v);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "\tTicks:" << std::setw(7) << float((end-start).count())
		<< " in " << name << std::endl;
}

template<int ROOT, typename T_Float, unsigned NewtonSteps>
void generate_root_functions(BEST_APPROX_BASIS basis)
{
	auto best = RootApprox_Best<ROOT, T_Float, NewtonSteps>(basis);
	
	std::cout << "/*" << std::endl;
	char name[] = "0 newtonian steps";
	name[0] = char('0' + NewtonSteps);
	Print_Test_Root_Approx(name, best, ROOT);
	std::cout << "*/" << std::endl;
	
	std::cout << best << std::endl << std::endl;
}

static float identity     (const float y)    {return y;}
static float std_sqrt     (const float y)    {return std::sqrt(y);}
static float std_sqrt_sqrt(const float y)    {return std::sqrt(std::sqrt(y));}
static float inverse      (const float y)    {return 1.f / y;}
static float inv_std_sqrt (const float y)    {return 1.f / std::sqrt(y);}
static float pow_quarter  (const float y)    {return std::pow(y,-.25f);}

int main(int argc, const char * argv[])
{
#if _MSC_VER
	// FTZ and DAZ flags for speed
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif


	std::cout << std::hex;
	
	for (auto &v : TEST_VALUES) v = (std::rand() / float(RAND_MAX));
	
	//auto std_sqrt = static_cast<float (*)(float)>(&std::sqrt);
	
	Print_Func_Profile("y",         identity);
	Print_Func_Profile("sqrt(x)", std_sqrt);
	Print_Func_Profile("sqrt(sqrt(y))", std_sqrt_sqrt);
	Print_Func_Profile("1/sqrt(y)", inv_std_sqrt);
	Print_Func_Profile("1/y",         inverse);
	Print_Func_Profile("pow(y,-.25)", pow_quarter);
	Print_Func_Profile("y",         identity);
	Print_Func_Profile("rb_2_root", rb_2_root);
	Print_Func_Profile("rb_3_root", rb_2_root);
	Print_Func_Profile("rb_4_root", rb_2_root);
	Print_Func_Profile("rb_inv_2_root", rb_inv_2_root);
	Print_Func_Profile("rb_inv_3_root", rb_inv_4_root);
	Print_Func_Profile("rb_inv_4_root", rb_inv_4_root);
	
	Print_Test_Root_Approx("std::sqrt", std_sqrt, 2);
	Print_Test_Root_Approx("rb_2_root",  rb_2_root,  2);
	Print_Test_Root_Approx("1/std::sqrt",  inv_std_sqrt,  -2);
	Print_Test_Root_Approx("rb_inv_2_root",  rb_inv_2_root,  -2);
	Print_Test_Root_Approx("rb_inv_4_root",  rb_inv_4_root,  -4);
	
	std::cout << std::endl << std::endl;

	std::cout << "#pragma once" << std::endl;
	std::cout << "#include <stdint.h>" << std::endl;
	std::cout << std::endl << std::endl;
	
	std::cout << "// Functions optimized for worst-case error" << std::endl;
	std::cout << std::endl << std::endl;
	
	generate_root_functions< 2,float,1>(BEST_WORST_CASE);
	generate_root_functions<-2,float,1>(BEST_WORST_CASE);
	generate_root_functions< 3,float,1>(BEST_WORST_CASE);
	generate_root_functions<-3,float,1>(BEST_WORST_CASE);
	generate_root_functions< 4,float,1>(BEST_WORST_CASE);
	generate_root_functions<-4,float,1>(BEST_WORST_CASE);
	
	/*generate_root_functions< 2,double,1>(BEST_WORST_CASE);
	generate_root_functions<-2,double,1>(BEST_WORST_CASE);
	generate_root_functions< 3,double,1>(BEST_WORST_CASE);
	generate_root_functions<-3,double,1>(BEST_WORST_CASE);
	generate_root_functions< 4,double,1>(BEST_WORST_CASE);
	generate_root_functions<-4,double,1>(BEST_WORST_CASE);*/
	
	//std::cout << "RootApprox<-2,float,1> error: " << std::flush;
	//std::cout << Test_RMS_Error(classicinvsqrt, -.5f, 1.f, 2.f) << std::endl;
	
	//std::cout << classicinvsqrt << std::endl;

	std::cin.ignore(255, '\n');
	
	return 0;
}
