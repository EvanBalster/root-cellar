#pragma once


#include <cstdint>
#include <cmath>
#include <algorithm>
#include <utility>
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
		struct pow_i_
		{
			template<typename X>
			static X calc(const X x)
			{
				if (EXP_INDEX > 0) return pow_i_<std::max(EXP_INDEX-1,0)>::calc(x) * x;
				if (EXP_INDEX < 0) return pow_i_<std::min(EXP_INDEX+1,0)>::calc(x) / x;
				return X(1);
			}
		};
		template<> struct pow_i_<-4> {template<typename X> static X calc(const X x) {return X(1)/((x*x)*(x*x));}};
		template<> struct pow_i_<-3> {template<typename X> static X calc(const X x) {return X(1)/(x*x*x);}};
		template<> struct pow_i_<-2> {template<typename X> static X calc(const X x) {return X(1)/(x*x);}};
		template<> struct pow_i_<-1> {template<typename X> static X calc(const X x) {return X(1)/x;}};
		template<> struct pow_i_< 0> {template<typename X> static X calc(const X x) {return X(1);}};
		template<> struct pow_i_< 1> {template<typename X> static X calc(const X x) {return x;}};
		template<> struct pow_i_< 2> {template<typename X> static X calc(const X x) {return x*x;}};
		template<> struct pow_i_< 3> {template<typename X> static X calc(const X x) {return x*x*x;}};
		template<> struct pow_i_< 4> {template<typename X> static X calc(const X x) {return (x*x)*(x*x);}};
		
		template<int ROOT_INDEX>
		struct root_i_
		{
			static_assert(ROOT_INDEX != 0, "0th root is invalid!");
			
			template<typename X>
			static X calc(const X x) {return std::pow(x, X(1)/X(ROOT_INDEX));}
		};
		template<> struct root_i_<-4> {template<typename X> static X calc(const X x) {return X(1)/std::sqrt(std::sqrt(x));}};
		template<> struct root_i_<-2> {template<typename X> static X calc(const X x) {return X(1)/std::sqrt(x);}};
		template<> struct root_i_<-1> {template<typename X> static X calc(const X x) {return X(1)/x;}};
		template<> struct root_i_< 1> {template<typename X> static X calc(const X x) {return x;}};
		template<> struct root_i_< 2> {template<typename X> static X calc(const X x) {return std::sqrt(x);}};
		template<> struct root_i_< 4> {template<typename X> static X calc(const X x) {return std::sqrt(std::sqrt(x));}};
	}
	
	template<int EXP_INDEX, typename T_Num>
	inline T_Num pow_i(const T_Num n)    {return detail::pow_i_<EXP_INDEX>::calc(n);}
	
	template<int ROOT_INDEX, typename T_Num>
	inline T_Num root_i(const T_Num n)    {return detail::root_i_<ROOT_INDEX>::calc(n);}
	
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
		double min_error = 0.0, min_error_arg = 0.0;
		double max_error = 0.0, max_error_arg = 0.0;
		
		double worst_error() const    {return std::max(-min_error, max_error);}
	};
	
	template<int ROOT_INDEX, typename T_Approx, typename T_Float>
	inline PowApprox_Stats Test_Root_Approx(
		const T_Approx &approx,
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
		measure_t sum_error = 0.0, sum_sq_error = 0.0,
			min_error     = 1e20, max_error     = -1e20,
			min_error_arg = 0.0, max_error_arg = 0.0;
		for (int i = ib; i <= ie; ++i)
		{
			float_t y = reinterpret_int_float(i), x = root_i<ROOT_INDEX>(y);
			measure_t error = (approx(y) - x) / x;
			sum_error += error;
			sum_sq_error += error*error;
			if (error < min_error) {min_error = error; min_error_arg = y;}
			if (error > max_error) {max_error = error; max_error_arg = y;}
		}
		double samples = double(ie - ib);
		return {
			sum_sq_error / samples,
			sum_error / samples,
			min_error, min_error_arg,
			max_error, max_error_arg};
	}
	
	/*template<typename T_Approx, typename T_Float>
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
	}*/
	
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
		//float_t exponent = float_t(1)/ROOT;
			
		union {float_t xf; int_t xi;};
		union {float_t yf; int_t yi;};
		/*yf = std::pow(range_min, exponent); ((ROOT > 0) ? --yi : ++yi);
		float_t xf_lower = detail::pow_i<ROOT>::calc(yf);*/
			
		// Measurements...
		float_t worst_error = 0.0;
		for (xi = ib; xi <= ie; ++xi)
		{
			// Move root-value to closest approximation
			yf = root_i<ROOT>(xf);
			/*while (true)
			{
				((ROOT > 0) ? ++yi : --yi);
				float_t xf_upper = detail::pow_i<ROOT>::calc(yf);
				if (xf_upper > xf)
				{
					if (xf_upper - xf > xf - xf_lower) ((ROOT > 0) ? --yi : ++yi);
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
		static const T_Real k = T_Real(N-1) / T_Real(N), p = T_Real(1)/T_Real(N);
		if (N > 0) return x *  k + p * y / pow_i<N-1>(x);
		else       return x * (k + p * y * pow_i<-N>(x));
	}
	
	/*
		A formula for a approximate roots affording fast implementation.
	 */
	template<int N, typename T_Float, unsigned NewtonSteps = 1>
	struct RootApprox
	{
		static_assert(N != 0, "0th root is invalid");
		
		static const int DEG = ((N>0) ? N : -N);
		
		using float_t  = T_Float;
		using range_t  = std::pair<float_t, float_t>;
		using as_int_t = float_as_int_t<float_t>;
		
		as_int_t constant;
		float_t  newton_m = float_t(1) / float_t(N);
		
		RootApprox(as_int_t _constant) :
			constant(_constant) {}
		
		/*
			Initial estimate based on float-reinterpretation hack.
		*/
		//static const as_int_t _rshift = as_int_t(std::ceil(std::log2(std::abs(N))));
		
		float_t initialEstimate(const float_t y) const
		{
			// Floating-point hack for initial estimate
			union {float_t x; as_int_t i;}; x = y;
			i = constant + i / as_int_t(N);
			/*if (DEG&(DEG-1)) i = constant + i / as_int_t(N);
			else if (N > 0)  i = constant + (i >> _rshift);
			else             i = constant - (i >> _rshift);*/
			return x;
		}
		float_t initialEstimate_inverse(const float_t x) const
		{
			// Approximate exponential; the inverse of the initial estimator
			union {float_t y; as_int_t i;}; y = x;
			i = (i - constant) * as_int_t(N);
			/*if (DEG&(DEG-1)) i = -constant + i * as_int_t(N);
			else if (N > 0)  i = -constant + (i << _rshift);
			else             i = -constant - (i << _rshift);*/
			return y;
		}
		/*
			One step of newtonian refinement.
		*/
		float_t newtonianRefinement(const float y, const float x) const
		{
			if (N > 0) return x *  (float_t(1)-newton_m) + newton_m * y / pow_i<N-1>(x);
			else       return x * ((float_t(1)-newton_m) + newton_m * y * pow_i<-N>(x));
		}
		
		/*
			Complete calculation.
		*/
		float_t operator()(const float_t y) const
		{
			float_t x = initialEstimate(y);
			
			for (unsigned i = 0; i < NewtonSteps; ++i)
				x = newtonianRefinement(y, x);
			
			return x;
		}
		
		// Get the suggested testing range for this root
		static range_t test_param_range()
		{
			return std::make_pair(float_t(1), float_t(1<<std::abs(N)));
		}
		
		/*
			Calculate range of relative error
		*/
		range_t errorRange_initial() const
		{
			const double P = double(1)/double(N);
			range_t range(float_t(1e20), float_t(-1e20));
			
			auto consider = [&](const float_t y)
			{
				float_t x = initialEstimate(y);
				float_t ratio = x / root_i<N>(y);
				range.first  = std::min(range.first,  ratio);
				range.second = std::max(range.second, ratio);
				//std::cout << "\tf(" << y << ") = " << x << " (/x = " << ratio << " - 1 = " << (ratio-1) << ")" << std::endl;
				return x;
			};
			
			// Locate the output-value discontinuity
			double ys = initialEstimate_inverse(float_t(1));
			if (ys < 1.0) ys *= double(1 << DEG);
			//std::cout << "\tys = " << ys << std::endl;
			
			// Interate over discontinuities
			double y1 = float_t(1), x1 = initialEstimate(y1);
			int i = 1;
			while (true)
			{
				// Select the next input or output discontinuity.
				double y2 = float_t(1<<i);
				if (y1 < ys && y2 > ys) y2 = ys;
				else                    ++i;
				
				// Consider the discontinuity itself
				double x2 = initialEstimate(y2);
				
				// Consider any local minimum or maximum
				{
					// TODO improve numerical stability here
					//float_t b = (x2-x1)/(y2-y1), a = x1 - b*y1;
					//float_t yM = -P*a / ((P-1)*b);
					
					// TODO this formula is wrong somehow
					float_t yM = (P/(P-1.0)) * (y1 - x1 * (y2-y1)/(x2-x1));
					//std::cout << "\t\tyM [" << y1 << "," << y2 << "] = " << yM << std::endl;
					if (yM > y1 && yM < y2)
					{
						//std::cout << 'M';
						consider(yM);
					}
				}
				
				//if (y2 == ys) std::cout << 'S';
				consider(y2);
				
				if (i > DEG) break;
				y1 = y2;
				x1 = x2;
			}
			
			return range;
		}
		range_t errorRange_refine(range_t prevRange) const
		{
			const float_t exponent = float_t(1)/float_t(N);
			range_t range(float_t(1e20), float_t(-1e20));
			
			auto consider = [&](const float_t ratio)
			{
				float_t refined
					= (1 - newton_m) * ratio
					+ newton_m * pow_i<1-N>(ratio);
				range.first  = std::min(range.first,  refined);
				range.second = std::max(range.second, refined);
				//std::cout << "Consider NR(" << ratio << ") = " << refined << std::endl;
			};
			consider(prevRange.first);
			consider(prevRange.second);
			
			// Consider additional local min/max.
			float_t extremum = root_i<N>(
				(newton_m * (exponent - float_t(1))) /
				(exponent * (newton_m - float_t(1))));
			if (extremum > prevRange.first && extremum < prevRange.second)
				consider(extremum);
			
			return range;
		}
		
		range_t errorRange() const
		{
			range_t range = errorRange_initial();
			
			for (unsigned i = 0; i < NewtonSteps; ++i)
				range = errorRange_refine(range);
				
			return range;
		}
		
		float_t error_worstCase() const
		{
			range_t range = errorRange();
			return std::max(std::abs(range.first-float_t(1)), std::abs(range.second-float_t(1)));
		}
	};
	
	template<typename I>
	I nextpow2(const I v)
	{
		I m = std::max<I>(v,1)-1;
		m = m | (v>>1) | (v>>2) | (v>>4);
		if (sizeof(I) > 1) m |= v >> 8;
		if (sizeof(I) > 2) m |= v >> 16;
		if (sizeof(I) > 4) m |= v >> 32;
		if (sizeof(I) > 8) m |= v >> 64;
		return m+1;
	}
	
	/*
	
	 */
	enum BEST_APPROX_BASIS
	{
		BEST_WORST_CASE  = 0,
		BEST_MEAN_SQUARE = 1,
		APPROX_WORST_CASE = 2,
	};
	 
	template<int N, typename T_Float, unsigned NewtonSteps = 1, BEST_APPROX_BASIS Basis = BEST_WORST_CASE>
	RootApprox<N,T_Float,NewtonSteps> RootApprox_Best()
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
			p           = float_t(1)/float_t(N),
			one_minus_p = float_t(1) - p;
		as_int_t
			k_min = as_int_t(std::floor(one_minus_p * L * (float_t(B) - sigma_max))),
			k_max = as_int_t(std::ceil (one_minus_p * L * (float_t(B) - sigma_min))),
			m_min = reinterpret_float_int(p),
			m_max = reinterpret_float_int(p*float_t(1.5));
		if (m_min > m_max) std::swap(m_min, m_max);
		if (NewtonSteps == 0) m_max = m_min;
			
		// Determine testing range...
		float_t
			test_min = float_t(1),
			test_max = float_t(1 << std::abs(N));
		
		auto get_score = [=](as_int_t k, as_int_t m) -> float_t
		{
			RootApprox<N, T_Float, NewtonSteps> candidate(k);
			candidate.newton_m = reinterpret_int_float(m);
			
			switch (Basis)
			{
			default:
			case BEST_WORST_CASE:   return std::abs(Test_Root_Approx_WorstCase<N>(candidate, test_min, test_max));
			case APPROX_WORST_CASE: return candidate.error_worstCase();
			case BEST_MEAN_SQUARE:  return float_t(Test_Root_Approx<N>(candidate, test_min, test_max).mean_sq_error);
			}
		};
		
		std::cout << std::hex << "//Searching k in [0x"
			<< k_min << ",0x" << k_max
			<< "], m in [" << reinterpret_int_float(m_min)
			<< "," << reinterpret_int_float(m_max) << "] ";
		
		float_t  best_score = float_t(1e20);
		as_int_t best_k = -1, best_m = -1;
		as_int_t
			k_lo = k_min, k_hi = k_max,
			m_lo = m_min, m_hi = m_max,
			k_step = (k_max - k_min) / 8,
			m_step = (m_max - m_min) / 8;
			
		k_step = nextpow2(k_step);
		m_step = nextpow2(m_step);
		k_step = m_step = std::max(k_step, m_step);
			
		while (k_lo < k_hi || m_lo < m_hi)
		{
			std::cout << '.' << std::flush;
			if (k_step == 0) k_step = 1;
			if (m_step == 0) m_step = 1;
			as_int_t
				k_start = k_lo + ((k_hi-k_lo)/k_step)/2,
				m_start = m_lo + ((m_hi-m_lo)/m_step)/2;
			for (as_int_t k = k_start; k <= k_hi; k += k_step)
				for (as_int_t m = m_start; m <= m_hi; m += m_step)
			{
				float_t score = get_score(k, m);
				
				if (score < best_score)
				{
					best_score = score;
					best_k = k;
					best_m = m;
				}
			}
			
			k_step = ((k_step > 1) ? std::max<as_int_t>(k_step>>2, 1) : 0);
			m_step = ((m_step > 1) ? std::max<as_int_t>(m_step>>2, 1) : 0);
			k_lo = std::max(k_min, best_k - 4 * k_step);
			k_hi = std::min(k_max, best_k + 4 * k_step);
			m_lo = std::max(m_min, best_m - 4 * m_step);
			m_hi = std::min(m_max, best_m + 4 * m_step);
		}
		
		std::cout << std::endl;
		
		
		/*as_int_t l = constant_min, r = constant_max;
		float_t l_score = get_score(l), r_score = get_score(r);
		while (l+1 < r)
		{
			// Score the value in the middle
			as_int_t m = l + (r-l)/2;
			float_t m_score = get_score(m);
			
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
			//std::cout << "//  ...exhaustive between " << l << " and " << r << "..." << std::endl;
			for (as_int_t k = l; k <= r; ++k)
			{
				float_t score = get_score(k);
				if (score < best_score)
				{
					best_score = score;
					best_constant = k;
				}
			}
		}*/
		RootApprox<N, T_Float, NewtonSteps> result(best_k);
		result.newton_m = reinterpret_int_float(best_m);
		
		std::cout << "//  ...best design k=" << best_k
			<< ", m=" << result.newton_m
			<< " with error score " << best_score << std::endl;
		return result;
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
		//static const float_t
		//	newton_m   = approx.newton_m, //float_t(N-1) / float_t(N),
		//	newton_div = float_t(1)   / float_t(absN);
		const int refine_power = absN - (N>0);
		for (unsigned i = 0; i < NewtonSteps; ++i)
		{
			out << "\tx" << ((N>0) ? " = " : " *= ")
				<< (float_t(1)-approx.newton_m) << float_suff << ((N>0) ? " * x" : "")
				<< ((N>0) ? " + " : " - ")
				<< ((N>0) ? approx.newton_m : -approx.newton_m) << float_suff << " * y";
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
