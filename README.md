
# Plotting Basins of a Univariate Rational map


This package PBURF allows us to visualize the attraction basins associated to the end points of a discrete semi-flow induced by a rational function on the Riemann sphere by using its geometry and complex structure.

The main advantage of the algorithm used by this package is that avoids the problem of overflows caused by denominators closed to zero and the problem of indetermination which appears when simultaneously the numerator and denominator are equal to zero. This is solved by working with homogeneous coordinates and the iteration of a homogeneous pair on the augmented complex projective line (Riemann spheres plus an additional superzero point).

This can be applied to any numerical method which construct a rational map to solve an univariate polynomial equation and verifying  that the roots of the equation are  fixed points of the associated rational map.


# Usage

To load the code just type

>using PBURF

The main function of package is

* `plottingBasinsUnivariateRationalFunctions(coefficientlistnum,coefficientlistden,...)`: This main function takes two lists of complex numbers that correspond to the coefficients of the numerator and denominator of a rational function and returns the list of fixed points of the function, a rectangular plot of the basins of these fixed points  and a color palette with colors associated to the different fixed points.

This function uses the following arguments

* `coefficientlistnum`: is the list of the complex coefficients of the numerator of a univariate rational function.

* `coefficientlistden`: is the list of the complex coefficients of the denominator of a univariate rational function                  .

This function has the following optional arguments

* `expresolution`: It is a nonnegative integer. The function gives a rectangular plot such that the sides are divided into 2^expresolution subintervals. It default value is expresolution=8.

* `iterationmax`: It is a nonnegative integer. This is the maximun of possible iterations of the rational function when you run this main function. It default value is iterationmax=25.

* `iterprecision`: It is a nonnegative integer. The stopping criterium is that chordal distance between two the last consecutive iterates is less than 1/10^iterationmax. It default value is iterationmax=3.

* `aproxprecision`: It is a nonnegative integer. When the chordal distance between the last iterate and a point of the list of fixed points is less than 1/10^aproxprecision the function assigns the starting point to the basin of this fixed point. It default value is aproxprecision=3.

This function also uses the following keyword arguments

* `colorstrategy`: It is an AbstractString. It posible values are "positionplusiteration", "iteration", and "positionfixedpoints". The defaut value is  colorstrategy="positionplusiteration". The program provides three different strategies for assigning colors to the different points of a plot.

* `model`:It is an AbstractString. The defaut value is  model="localrectangle". Future development of this program will use other models like the Riemann sphere, a pair of disks, etc.


* `rectanglesides`: It is a 4-dimensional tuple of real numbers (a,b,c,d)  which corresponds to rectangle obtained as the product of  the intervarls  [a,b] and [c,d]. The defaut values is  rectanglesides=(-1.5,1.5,-1.5,1.5).



#Installation

This software can be instaled  by giving the following command in the julia command line
(or in some [IJulia](https://github.com/JuliaLang/IJulia.jl) notebook):

> Pkg.clone("https://github.com/luisjavierhernandez/PBURF.jl")

in this way, all dependencies will be satisfied automatically.

The code will be upgraded every time the Pkg.update() command is used.

Alternatively, you can manually copy the whole directory structure
to your julia package directory (use Pkg.dir() to locate it),
and then run Pkg.update() to download the dependencies.


You will need to have the Python [Matplotlib](http://matplotlib.org/)
library installed on your machine in order to use PyPlot.  You can either
do inline plotting with [IJulia](https://github.com/JuliaLang/IJulia.jl),
which doesn't require a GUI backend, or use the Qt, wx, or GTK+ backends
of Matplotlib.

## Dependencies

*Polynomials.jl

*Pyplot.jl

*Colorsjl
