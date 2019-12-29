# Root Cellar

This library is an experiment in generalizing the famed "fast inverse square root" hack to other exponents.  For example:

* Fast cubic root
* Fast inverse 4th root
* More/less precise versions of fast roots
* Fast roots optimized for mean-square error instead of worst-case error

In time, this library will provide a reference table of optimal constants, precision and speed measurements for many different variations on the formulas.



## Notes

I'm researching fast roots for applications in signal processing and graphics rendering — and as a fun distraction from more intensive research work.

Fast roots are useful in performance-intensive applications where a small amount of relative error is tolerable.  Modern CPUs offer intrinsic support for `sqrt(y)`, making a fast approximation unnecessary.  The other approximations in this library often outperform their standard equivalents, however, including `1/sqrt(y)` and especially calls to `pow`.



## Anatomy of the Algorithms

The fast root functions studied with this library have the following steps:

* Reinterpret the bits of floating-point parameter `y` as an integer `i`.
* Set `i` to `k + i * p` where `p` is a fractional exponent and `k` is a "magic constant".
* Reinterpret the new `i` as a floating-point value `x`.
* Apply zero or more refinements to `x`
  * Refinements are based on the Newtonian formula:
    `x *= (1-k) + k * y / x^(1/p)`
    Where `k` is close to `p`.

This library supports only roots — exponents `p = 1/N` where `N` is some integer.



## Assessing Formulas

There is some flexibility in the design of fast approximate roots.  Historic implementations have been tuned to minimize worst-case error.  Applications exist, however, where it may be better to minimize mean error, mean-squared error or error when `y` is very close to 1.

This library specifies root approximations with four parameters:

* **The exponent `p`** *(1/N for integer N)*
* **The "magic constant" `k`** *(integer)*
* **Number of Newtonian refinement steps** *[0,∞)*
* **Newtonian refinement constant** *(floating-point)*

Designs are evaluated using these criteria:

* Relative computation time
* Root-mean-square relative error `sqrt(mean: (approx_f(y) - f(y))^2 / y)`
* Mean relative error `mean: (approx_f(y) - f(y)) / y`
* Worst-case relative error `(approx_f(y) - f(y)) / y`
* *Planned*: Relative error for values near 1



## References

[0x5f3759df (Christian Plesner Hansen)](http://h14s.p5r.org/2012/09/0x5f3759df.html) on the theory behind the hack, derivation of the "magic constant" and generalizing the approximation to exponents between -1 and 1.

[Understanding Quake's Fast Inverse Square Root (BetterExplained)](https://betterexplained.com/articles/understanding-quakes-fast-inverse-square-root) includes an explanation of the algebra behind the Newtonian refinement used in the hack.

[Inverse Square Root (University of Waterloo)](https://ece.uwaterloo.ca/~dwharder/aads/Algorithms/Inverse_square_root/) describes how we can reduce worst-case error by adjusting the multipliers used in Newton's method.



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