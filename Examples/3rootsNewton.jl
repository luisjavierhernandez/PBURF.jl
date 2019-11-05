
import Pkg

Pkg.status()
using PyCall
pygui(:tk)
using PyPlot

import PBURF

coefficientlistnumcompila=complex([1.,0.,0.,2.])
coefficientlistdencompila=complex([0.,0.,3.,0.])
result=PBURF.plottingBasinsUnivariateRationalFunctions(coefficientlistnumcompila,
    coefficientlistdencompila)

gcf()

PBURF.plottingBasinsUnivariateRationalFunctions(coefficientlistnumcompila,
    coefficientlistdencompila; colorstrategy="iteration")

gcf()
