#pragma once


#include <cstdint>
#include <cmath>
#include <ostream>

#include <iostream> //debug


namespace rootbeer
{
	namespace detail
	{
		template<typename T_Float> struct float_traits {};
		template<> struct float_traits<float>
		{
			using as_int_t = int32_t;
			static const as_int_t
				bits_exponent =  8,
				bits_mantissa = 23;
			static const char *name() {return "float";}
			static const char *suffix() {return "f";}};
		template<> struct float_traits<double>
		{
			using as_int_t = int64_t;
			static const as_int_t
				bits_exponent = 11,
				bits_mantissa = 52;
			static const char *name() {return "double";}
			static const char *suffix() {return "";}
		};
		
		template<typename T_Real> struct int_traits {};
		template<> struct int_traits<int32_t> {using as_real_t = float; static const char *name() {return "int32_t";}};
		template<> struct int_traits<int64_t> {using as_real_t = double; static const char *name() {return "int64_t";}};
		template<> struct int_traits<uint32_t> {using as_real_t = float; static const char *name() {return "uint32_t";}};
		template<> struct int_traits<uint64_t> {using as_real_t = double; static const char *name() {return "uint64_t";}};
		
		template<int EXP_INDEX>
		struct pow_i
		{
			template<typename X>
			static X calc(const X x)
			{
				if (EXP_INDEX > 0) return pow_i<std::max(EXP_INDEX-1,0)>::calc(x) * x;
				if (EXP_INDEX < 0) return pow_i<std::min(EXP_INDEX+1,0)>::calc(x) / x;
				return X(1);
			}
		};
		template<> struct pow_i<-4> {template<typename X> static X calc(const X x) {return X(1)/((x*x)*(x*x));}};
		template<> struct pow_i<-3> {template<typename X> static X calc(const X x) {return X(1)/(x*x*x);}};
		template<> struct pow_i<-2> {template<typename X> static X calc(const X x) {return X(1)/(x*x);}};
		template<> struct pow_i<-1> {template<typename X> static X calc(const X x) {return X(1)/x;}};
		template<> struct pow_i< 0> {template<typename X> static X calc(const X x) {return X(1);}};
		template<> struct pow_i< 1> {template<typename X> static X calc(const X x) {return x;}};
		template<> struct pow_i< 2> {template<typename X> static X calc(const X x) {return x*x;}};
		template<> struct pow_i< 3> {template<typename X> static X calc(const X x) {return x*x*x;}};
		template<> struct pow_i< 4> {template<typename X> static X calc(const X x) {return (x*x)*(x*x);}};
		
	}
	
	template<typename T_Float>
	using float_as_int_t = typename detail::float_traits<T_Float>::as_int_t;
	template<typename T_Int>
	using int_as_float_t = typename detail::int_traits<T_Int>::as_real_t;
	
	template<typename T_Float>
	float_as_int_t<T_Float> reinterpret_float_int(const T_Float v)    {return * reinterpret_cast<const float_as_int_t<T_Float>*>(&v);}
	template<typename T_Int>
	int_as_float_t<T_Int>   reinterpret_int_float(const T_Int   v)    {return * reinterpret_cast<const int_as_float_t<T_Int>*>(&v);}

	/*
		Calculate the root-mean-square error of an exponent approximation.
	 */
	struct PowApprox_Stats
	{
		double mean_sq_error  = 0.0;
		double mean_error = 0.0;
		double worst_error = 0.0;
	};
	
	template<typename T_Approx, typename T_Float>
	inline PowApprox_Stats Test_Pow_Approx(
		const T_Approx &approx,
		T_Float exponent,
		T_Float range_min,
		T_Float range_max)
	{
		using float_t = T_Float;
		using int_t = float_as_int_t<float_t>;
		int_t
			ib = reinterpret_float_int(range_min),
			ie = reinterpret_float_int(range_max);
			
		// Measurements...
		using measure_t = float_t;
		measure_t sum_error = 0.0, sum_sq_error = 0.0, worst_error = 0.0;
		for (int i = ib; i <= ie; ++i)
		{
			float_t x = reinterpret_int_float(i), y = std::pow(x, exponent);
			measure_t error = (approx(x) - y) / y;
			sum_error += error;
			sum_sq_error += error*error;
			if (std::abs(error) > std::abs(worst_error)) worst_error = error;
		}
		double samples = double(ie - ib);
		return {
			sum_sq_error / samples,
			sum_error / samples,
			worst_error};
	}
	
	template<int ROOT, typename T_Approx, typename T_Float>
	inline float Test_Root_Approx_WorstCase(
		const T_Approx &approx,
		T_Float range_min,
		T_Float range_max)
	{
		using float_t = T_Float;
		using int_t = float_as_int_t<float_t>;
		int_t
			ib = reinterpret_float_int(range_min),
			ie = reinterpret_float_int(range_max);
		float_t
			exponent = float_t(1)/ROOT;
			
		union {float_t xf; int_t xi;};
		union {float_t yf; int_t yi;};
		yf = std::pow(range_min, float_t(1)/ROOT);
		--yi;
		//float_t xf_lower = detail::pow_i<ROOT>::calc(yf);
			
		// Measurements...
		float_t worst_error = 0.0;
		for (xi = ib; xi <= ie; ++xi)
		{
			// Move root-value to closest approximation
			yf = std::pow(xf, exponent);
			/*while (true)
			{
				++yi;
				float_t xf_upper = detail::pow_i<ROOT>::calc(yf);
				if (xf_upper > xf)
				{
					if (xf_upper - xf > xf - xf_lower) --yi;
					break;
				}
				xf_lower = xf_upper;
			}*/
			
			float_t error = (approx(xf) - yf) / yf;
			//sum_error += error;
			//sum_sq_error += error*error;
			if (std::abs(error) > std::abs(worst_error)) worst_error = error;
		}
		return worst_error;
	}
	
	
	/*
		Newtonian step for refining x toward the Nth root of y
	 */
	template<int N, typename T_Real>
	T_Real newtonian_for_root(const T_Real x, const T_Real y)
	{
		static const T_Real k = T_Real(N-1) / T_Real(N), invN = T_Real(1)/T_Real(N);
		if (N > 0) return x *  k + invN * y / detail::pow_i<N-1>::calc(x);
		else       return x * (k + invN * y * detail::pow_i<-N>::calc(x));
	}
	
	/*
		A formula for a approximate roots affording fast implementation.
	 */
	template<int N, typename T_Float, unsigned NewtonSteps = 1>
	struct RootApprox
	{
		static_assert(N != 0, "0th root is invalid");
		
		using float_t  = T_Float;
		using as_int_t = float_as_int_t<float_t>;
		
		as_int_t constant;
		const as_int_t _rshift;
		
		RootApprox(as_int_t _constant) :
			constant(_constant),
			_rshift(std::ceil(std::log2(std::abs(N)))) {}
		
		float_t operator()(const float_t y) const
		{
			// Floating-point hack for initial estimate
			union {float_t x; as_int_t i;};
			x = y;
			if (std::abs(N)&(std::abs(N)-1)) i = constant + i / as_int_t(N);
			else if (N > 0)                  i = constant + (i >> _rshift);
			else                             i = constant - (i >> _rshift);
			
			// Newtonian refinements
			for (unsigned i = 0; i < NewtonSteps; ++i)
				x = newtonian_for_root<N>(x, y);
			
			return x;
		}
	};
	
	/*
	
	 */
	enum BEST_APPROX_BASIS
	{
		BEST_WORST_CASE = 0,
		BEST_MEAN_SQUARE = 1,
	};
	 
	template<int N, typename T_Float, unsigned NewtonSteps = 1>
	RootApprox<N,T_Float,NewtonSteps> RootApprox_Best(BEST_APPROX_BASIS basis = BEST_WORST_CASE)
	{
		using float_t = T_Float;
		using as_int_t = float_as_int_t<float_t>;
		
		// Worst-case value of x - log2(1 + x)
		float_t
			sigma_min = float_t(.00000),
			sigma_max = float_t(.08608);
		
		as_int_t
			L = (as_int_t(1) << detail::float_traits<float_t>::bits_mantissa),
			B = (as_int_t(1) << (detail::float_traits<float_t>::bits_exponent-1)) - 1;
		float_t
			p           = float_t(1)/N,
			one_minus_p = float_t(1) - p;
		as_int_t
			constant_min = std::floor(one_minus_p * L * (float_t(B) - sigma_max)),
			constant_max = std::ceil (one_minus_p * L * (float_t(B) - sigma_min));
			
		// Determine testing range...
		float_t
			test_min = float_t(1),
			test_max = float_t(1 << std::abs(N));
		
		double   best_score = 1e40;
		as_int_t best_constant = -1;
		
		//std::cout << "Searching constants between 0x"
		//	<< constant_min << " and 0x" << constant_max << "..." << std::endl;
		
		auto get_score = [=](as_int_t k)
		{
			RootApprox<N, T_Float, NewtonSteps> candidate(k);
			switch (basis)
			{
			default:
			case BEST_WORST_CASE:  return std::abs(Test_Root_Approx_WorstCase<N>(candidate, test_min, test_max));
			case BEST_MEAN_SQUARE: return float(Test_Pow_Approx(candidate, p, test_min, test_max).mean_sq_error);
			}
		};
		
		as_int_t l = constant_min, r = constant_max;
		double l_score = get_score(l), r_score = get_score(r);
		while (l+1 < r)
		{
			// Score the value in the middle
			as_int_t m = l + (r-l)/2;
			double m_score = get_score(m);
			
			if (m_score < best_score)
			{
				best_score = m_score;
				best_constant = m;
			}
			
			if (m_score > l_score && m_score > r_score) break;
			
			if (l_score < r_score)
			{
				r = m;
				r_score = m_score;
			}
			else
			{
				l = m;
				l_score = m_score;
			}
		}
		
		// Search a small area exhaustively
		l -= 256;
		r += 256;
		
		if (l+1 < r)
		{
			//std::cout << "    ...exhaustive between " << l << " and " << r << "..." << std::endl;
			for (as_int_t k = l; k <= r; ++k)
			{
				double score = get_score(k);
				if (score < best_score)
				{
					best_score = score;
					best_constant = k;
				}
			}
		}
		
		//std::cout << "    ...best constant is " << best_constant << " with error score " << best_score << std::endl;
		return RootApprox<N, T_Float, NewtonSteps>(best_constant);
	}
	
	
	
	
	/*
		Generalized newtonian step for Nth root:
		.    error(x)     x^N - y
		x =  --------  =  ---------
		.    error'(x)    Nx^(N-1)
	 */
}

template<int N, typename T_Float, unsigned NewtonSteps>
std::ostream &operator<<(std::ostream &out,
	const rootbeer::RootApprox<N, T_Float, NewtonSteps> &approx)
{
	static_assert(N != 0, "0th root is invalid");
	
	static const int absN = ((N<0)?-N:N);

	using float_t = T_Float;
	const char *float_decl = rootbeer::detail::float_traits<float_t>::name();
	const char *float_suff = rootbeer::detail::float_traits<float_t>::suffix();
	using as_int_t = rootbeer::float_as_int_t<float_t>;
	const char *int_decl = rootbeer::detail::int_traits<as_int_t>::name();
	
	out << std::hex;
	out << float_decl << " rb_";
	if (N < 0) out << "inv_";
	out << absN << "_root(const " << float_decl << " y)\n";
	out << "{\n";
	out << "\tunion {" << float_decl << " x; " << int_decl << " i;}; x = y; // interpret float as integer\n";
	
	// Magic line
	out << "\ti = 0x" << approx.constant << ((N>0) ? " + " : " - ")
		<< "(i";
	if (absN & (absN-1)) out << " / " << absN;
	else                 out << " >> " << int(std::log2(absN));
	out << "); // log-approximation hack\n";
	
	// Newtonian lines
	if (NewtonSteps)
	{
		static const float_t
		newton_k   = float_t(N-1) / float_t(N),
		newton_div = float_t(1)   / float_t(absN);
		const int refine_power = absN - (N>0);
		for (unsigned i = 0; i < NewtonSteps; ++i)
		{
			out << "\tx" << ((N>0) ? " = " : " *= ")
				<< newton_k << float_suff << ((N>0) ? " * x" : "")
				<< ((N>0) ? " + " : " - ")
				<< newton_div << float_suff << " * y";
			if (refine_power != 0)
			{
				out << ((N>0) ? " / " : " * ");
				
				if (refine_power == 1) out << "x";
				else
				{
					out << "(x";
					for (int i = 1; i < refine_power; ++i) out << "*x";
					out << ")";
				}
				out << "; // newtonian step #" << float(i+1) << "\n";
			}
		}
	}
	out << "\treturn x;\n";
	
	out << "}";
	
	return out;
}
