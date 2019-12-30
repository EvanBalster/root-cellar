# Root Cellar

This library is an experiment in generalizing the famed "fast inverse square root" hack to other exponents.  For example:

* Fast cubic root
* Fast inverse 4th root
* More/less precise variations
* Fast roots optimized for mean-square error instead of worst-case error

Fast roots are useful in performance-intensive applications where a small amount of relative error is tolerable.  Modern CPUs offer intrinsic support for `sqrt(y)`, making a fast approximation unnecessary.  The other approximations in this library often outperform their standard equivalents, however, including `1/sqrt(y)` and especially calls to `pow`.



## The Algorithm

This library defines root approximations in terms of four parameters:

- **Root Index** `N` *(integer)*
- **Magic Constant `K`** *(integer)*
- **Number of Newtonian refinement steps** `R` *(unsigned integer)*
- **Pseudo-Newtonian Constant** `M` *(floating-point, close to 1/N)*

The Nth root of a floating-point value `y` is approximated using these steps:

1. Reinterpret the bits of `y` as an integer `i`.
2. Calculate `K + i / N` and reinterpret it into a floating-point `x`, our initial estimate.
3. Apply `R` pseudo-Newtonian refinements, improving our estimate's accuracy.
   *  `x *= (1-M) + M * y / x^(1/p)`
4. Return `x`.

For example:

```c++
float fastinvsqrt(const float y)
{
/* 1 */  union {float x; int32_t i;}; x = y;
/* 2 */  i = 0x5f32a121 - (i >> 1);
/* 3 */  x *= 1.5351f - 0.535102f * y * (x*x);
/* 4 */  return x;
}
```

For more information about how this works, see the **References** section below.



## Assessing Formulas

Fast roots are conventionally optimized to minimize worst-case relative error.  Applications exist, however, where it may be better to minimize mean error, mean-squared error or error when `y` is very close to 1.

In this library, designs are evaluated using these criteria:

- Relative computation time
- Root-mean-square relative error `sqrt(mean: (approx_f(y) - f(y))^2 / y)`
- Mean relative error `mean: (approx_f(y) - f(y)) / y`
- Worst-case relative error `(approx_f(y) - f(y)) / y`



## Table of Constants

This table lists the best constants I've found for Nth roots with a single pseudo-Newtonian refinement.  For example, `N=2` corresponds to the square root while `N=-2` corresponds to the inverse square root.

All error values in this table are based on the relative error `(f(y) - y^(1/N)) / y^(1/N)`.  The error for which the function was optimized is marked in bold.

#### 32-bit floating-point roots

Error values were calculated exhaustively for these constants; I believe these are the most accurate designs for this parameterization.

| N    | 0 Refinements                                                | 1 Refinement                                                 | 2 Refinements                                                |
| ---- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| +2   | k: `0x1fbb4f2e`<br/>m: *n/a*<br/>eMax: **.0347475**<br/>eRMS: .0190506<br/>eMean: –.005133 | k: `0x1fbed49a`<br>m: `+.510929`<br/>
eMax: **.000239058**<br/>
eRMS: .000148007<br/>
eMean: –.0000454788 | k: `0x1fbb75ad`<br>m: `+.500122`<br/>
eMax: **1.68567e-7**<br/>
eRMS: 4.7967e-8<br/>
eMean: –1.41369e-8 |
| –2   | k: `0x5f37642f`<br/>m: *n/a*<br/>
eMax: **.0342129**<br/>eRMS: .0244769<br/>
eMean: +0.01338 | k: `0x5f32a121`<br/>m: `-.535102`<br/>
eMax: **.000773445**<br/>
eRMS: .000494072<br/>
eMean: –.000025526 | k: `0x5f3634f9`<br/>m: `-.501326`<br/>
eMax: **1.40452e-6**<br/>
eRMS: 8.95917e-7<br/>
eMean: +6.21959e-7 |
| +3   | k: `0x2a510680`<br/>m: *n/a* <br/>
eMax: **.0315547**<br/>
eRMS: .0180422<br/>
eMean: +.00321401 | k: `0x2a543aa3`<br/>m: `+.347252`<br/>
eMax: **.000430098**<br/>
eRMS: .000237859<br/>
eMean: –.0000977778 | k: `0x`<br/>m: `+.333818`<br/>
eMax: 6.45394e-7<br/>
eRMS: 2.90881e-7<br/>
eMean: –2.40255e-7 |
| –3   | k: `0x54a232a3`<br/>m: *n/a*<br/>
eMax: **.0342405**<br/>
eRMS: .0195931<br/>
eMean: +.0068072 | k: `0x549da7bf`<br/>m: `-.364707`<br/>
eMax: **.00102717**<br/>
eRMS: .000742809<br/>
eMean: +.000277928 | k: `0x54a1b99d`<br/>m: `-.334677`<br/>
eMax: **2.18458e-6**<br/>
eRMS: 1.04454e-6<br/>
eMean: +6.17437e-7 |
| +4   | k: `0x2f9b374e`<br/>m: *n/a*<br/>
eMax: **.0342323**<br/>
eRMS: .015625<br/>
eMean: +.00425431 | k: `0x2f9ed7c0`<br/>m: `+.266598`<br/>
eMax: **.000714053**<br/>
eRMS: .000444122<br/>
eMean: –.000195135 | k: `0x2f9b8068`<br/>m: `+.250534`<br/>
eMax: **9.49041e-7**<br/>
eRMS: 5.28477e-7<br/>
eMean: –4.76837e-07 |
| –4   | k: `0x4f58605b`<br/>m: *n/a*<br/>
eMax: **.0312108**<br/>
eRMS: .020528<br/>
eMean: +.00879536 | k: `0x4f542107`<br/>m: `-.277446`<br/>
eMax: **.00110848**<br/>
eRMS: .000733642<br/>
eMean: +.000175304 | k: `0x4f58020d`<br/>m: `-.251282`<br/>
eMax: **2.76944e-6**<br/>
eRMS: 1.3487e-6<br/>
eMean: +5.97779e-7 |

See `root_cellar_generated.h` for reference implementations of these functions.

#### 64-bit floating-point roots

These designs were found using a fast approximation of the maximum error.  There is room to improve these constants (with extremely marginal benefit).

| N    | 0 Refinements                                               | 1 Refinement                                                 | 2 Refinements                                                |
| ---- | ----------------------------------------------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| +2   | k: `0x1ff769e5b00cb024`<br/>m: *n/a*<br/>eMax: **.0347474** | k: `0x1ff7da9258189b10`<br />m: `+.51093`<br/>eMax: **.000238945** | k: `0x1ff76e33f8e94831`<br />m: `+.500124`<br/>eMax: **3.08405e-8** |
| –2   | k: `0x5fe6ec85e7de30da`<br/>m: *n/a*<br/>eMax: **.0342128** | k: `0x5fe65423e81eece9`<br/>m: `-.535103`<br/>
eMax: **.00077328** | k: `0x5fe6bbf0c11e182d`<br/>m: `-.501434`<br/>
eMax: **1.36764e-6** |
| +3   | k: `0x2a9f76253119d328`<br/>m: *n/a*<br/>
eMax: **.0315546** | k: `0x2a9fdca8d39b1833`<br/>m: `+.347251`<br/>
eMax: **.000429969** | k: `0x2a9f5317d3f76c27`<br/>m: `+.333791`<br/>
eMax: **4.7027e-7** |
| –3   | k: `0x553ef0ff289dd794`<br/>m: *n/a*<br/>
eMax: **.0342405** | k: `0x553e5fa2bf4bb94e`<br/>m: `-.364707`<br/>
eMax: **.001027** | k: `0x553eb1a359e5ec49`<br/>m: `-.335169`<br/>
eMax: **3.77555e-6** |
| +4   | k: `0x2ff366e9846f3cf9`<br/>m: *n/a*<br/>
eMax: **.0342321** | k: `0x2ff3daf850a16998`<br/>m: `+.266598`<br/>
eMax: **.00071393** | k: `0x2ff3578de1c1dc42`<br/>m: `+.250729`<br/>
eMax: **1.41358e-6** |
| –4   | k: `0x4feb0c0b7fa996ad`<br/>m: *n/a*<br/>
eMax: **.0312107** | k: `0x4fea8420dfe0c1b2`<br/>m: `-.277446`<br/>
eMax: **.0011083** | k: `0x4feaff5406bb3437`<br/>m: `-.251281`<br/>
eMax: **2.61417e-6** |



## Further Notes

I decided to research fast roots for applications in signal processing and graphics rendering — and as a fun distraction from more intensive research work.  I got in *way* over my head.

It is possible to use methods like the ones presented here to estimate any fixed exponent, such as `y^(3/2)` or `y^8`.  It is also possible to create an approximate `pow` function with a similar mechanism, and one can be found in the fastapprox library.

When using more than one Newtonian refinement, more accurate designs may be possible by varying the pseudo-Newtonian constant used for each refinement.  I have not explored this possibility as 3-dimensional configurations or higher would require more sophisticated solvers or a lot of computation time.



## References

[0x5f3759df (Christian Plesner Hansen)](http://h14s.p5r.org/2012/09/0x5f3759df.html) on the theory behind the hack, derivation of the "magic constant" and generalizing the approximation to exponents between -1 and 1.

[Understanding Quake's Fast Inverse Square Root (BetterExplained)](https://betterexplained.com/articles/understanding-quakes-fast-inverse-square-root) includes an explanation of the algebra behind the Newtonian refinement used in the hack.

[Inverse Square Root (University of Waterloo)](https://ece.uwaterloo.ca/~dwharder/aads/Algorithms/Inverse_square_root/) describes how we can reduce worst-case error by using a modified Newtonian refinement with a baked-in multiplier.

[Improving the Accuracy of the Fast Inverse Square Root Algorithm (Walczyk et al)](<https://www.researchgate.net/publication/323276298_Improving_the_accuracy_of_the_fast_inverse_square_root_algorithm>) discusses additional techniques for improving accuracy, beyond the ones presented here.



## License

```
Copyright 2019 Evan Balster

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```