# Root Cellar: Methodology

In these notes we refer to the parameters of a fast approximate root function as `y` and the result as `x`.  The fast root function approximates an exponential function `x = y^p` where `p = 1/N` and `N` is the root index.  For example, fast inverse square root would specify `N = -2, p = -.5`.



## Finding the Maximum Relative Error

For any exponent `p`, the integer-based approximation used in this library's past roots exhibits a relative error which is periodic with respect to the logarithm of `y`.

The integer-based approximation itself can be seen to be a piece-wise linear approximation of the root — that is, a sectioned function where each section can be formulated as `x = a + b*y`  Minimums and maximums in the relative error may occur:

* At the borders of each section (wherever `x` or `y` equals a power of two);
* At most once within each section, where `y = -pa/((p-1)b)`.

The Newtonian refinement formula is expressed in terms of exponent `p` and tweak constant `m`, where `m=p` for the canonical Newtonian formula.

It may be demonstrated to behave as a function of the previous relative error.  This formula is most simply expressible in terms of the error-ratio `r = x / y^p`:

> `r_new = r * (1 - m) + m * r / r^(1/p)`

The range of relative error may be determined by examining three cases:

* The minimum and maximum values of `r` from the previous step.
* The local maximum `r = (m*(p-1))/(p*(m-1))^p`, if between the previous minimum and maximum.

Combining our knowledge about minima and maxima in these steps, we can quickly evaluate the relative error of any approximate root or fixed-power function.  This quick evaluation allows us to quickly search for optimal parameters.