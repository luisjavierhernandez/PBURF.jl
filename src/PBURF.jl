module PBURF

import Polynomials
using Polynomials
import PyPlot;
#using PyPlot
import Colors
#using Colors

function homogeneousNormalization(
        twotuple::Tuple{Complex{Float64},Complex{Float64}},
        aproxprecision::Int64=15)
    tt1=complex(twotuple[1])
    tt2=complex(twotuple[2])
    if (abs(tt1)<1.0/10^aproxprecision && abs(tt2)<1.0/10^aproxprecision)
        hpoint=(complex(0.0),complex(0.0))
    else
        if abs(tt1)<= abs(tt2)
        hpoint=(tt1/tt2,complex(1.0))
        else
        hpoint=(complex(1.0),tt2/tt1)
        end
    end
    return hpoint
end

#Example:
homogeneousNormalization((1.0+0.0*im, 0.0+1.0*im))

function sphereBijection(twotuple::Tuple{Complex{Float64},Complex{Float64}},
                        aproxprecision::Int64=8)
    z=twotuple[1]
    t=twotuple[2]
    if (abs(z)<1.0/10^aproxprecision && abs(t)<1.0/10^aproxprecision)
    point=[0.0,0.0,0.0]
    else
    point=[real((conj(z)*t + conj(t)*z)/(conj(t)*t + conj(z)*z)),
            real((1im*(conj(z)*t - conj(t)*z))/(conj(t)*t + conj(z)*z)),
            real((-conj(t)*t + conj(z)*z)/(conj(t)*t + conj(z)*z))]
    end
    return point
end

sphereBijection((1.0+0.0*im, 0.0+1.0*im),9)

function chordalMetric(twotuple::Tuple{Complex{Float64},Complex{Float64}},
                        twotuple1::Tuple{Complex{Float64},Complex{Float64}},
                        aproxprecision::Int64=8)
    norma=vector->sqrt(vector[1]^2+vector[2]^2+vector[3]^2)
    return norma((sphereBijection(twotuple,aproxprecision)-
    sphereBijection(twotuple1,aproxprecision)))
end

chordalMetric((1.0+0.0*im, 0.0+1.0*im),(0.0+0.0*im, 1.0+0.0*im),9)

function bivariatepolyfunction(coefficientlist::Array{Complex{Float64},1},
     u::Complex{Float64}, t::Complex{Float64},d::Int64)
    ff=coefficientlist[1]*t^d
    for i in 2:length(coefficientlist)
        ff=ff+coefficientlist[i]*u^(i-1)*t^(d-i+1)
        +i
    end
    return ff
end


function bivariatepolyfunction(coefficientlist::Array{T,1},
    u::T, t::T,d::Int64) where {T<:Number}
    ff=coefficientlist[1]*t^d
    for i in 2:length(coefficientlist)
        ff=ff+coefficientlist[i]*u^(i-1)*t^(d-i+1)
        +i
    end
    return complex(ff)
end


function bivariatepolyfunction(coefficientlist::Array{T,1},
         u::Complex{Float64},
        t::Complex{Float64},d::Int64) where {T<:Number}
    coefficientlistcomplex=complex(coefficientlist)
    ff=coefficientlistcomplex[1]*t^d
    for i in 2:length(coefficientlist)
        ff=ff+coefficientlist[i]*u^(i-1)*t^(d-i+1)
    +i
    end
    return ff
end



function pairoffunctions(coefficientlistnum::Array{T,1},
    coefficientlistden::Array{T,1}) where {T<:Number}
    ln=length(coefficientlistnum)-1
    ld=length(coefficientlistden)-1
    d=max(ln,ld)
    fff(u::Complex{Float64}, t::Complex{Float64})=
    bivariatepolyfunction(coefficientlistnum,u,t,d)
    ggg(u::Complex{Float64}, t::Complex{Float64})=
    bivariatepolyfunction(coefficientlistden,u,t,d)
    return fff, ggg
end


coefficientlistnumcompila=complex([1.,0.,0.,2.])
coefficientlistdencompila=complex([0.,0.,3.,0.])

bivarej=bivariatepolyfunction(coefficientlistnumcompila,
    1.0+im* 2.3, 2.0+im*0.9,8)

hpaircompila=pairoffunctions(coefficientlistnumcompila,
    coefficientlistdencompila)

function rationalFunction(hpair::Tuple{Function,Function},
    twotuple::Tuple{Complex{Float64},Complex{Float64}},
    aproxprecision::Int64=8)
    c=twotuple[1]
    d=twotuple[2]
    F=hpair[1]
    G=hpair[2]
    cnew=F(c,d)
    dnew=G(c,d)
    hresult=homogeneousNormalization((cnew,dnew),aproxprecision)
    return  hresult
end

rationalFunction(hpaircompila,(1.0+0.0*im, 0.0+1.0*im),9)

function fixedPointsofaIrreduciblePair(
        coefficientlistnum::Union{Vector{Float64},Vector{Complex{Float64}}},
        coefficientlistden::Union{Vector{Float64},Vector{Complex{Float64}}},
        aproxprecision::Int64=8)
    coefficientlistnumcomplex=complex(coefficientlistnum)
    coefficientlistdencomplex=complex(coefficientlistden)
    if length(coefficientlistnum)!=length(coefficientlistden)
    println("The length of the first list
            is different of the length of the second")
    elseif  length(coefficientlistnum)==1 && length(coefficientlistden)==1
    finallistoffixedpoint=homogeneousNormalization((coefficientlistnumcomplex[1],
    coefficientlistdencomplex[1]),
    aproxprecision)
    return [finallistoffixedpoint]
    elseif abs(coefficientlistnum[1])==0.0 &&
        coefficientlistden[length(coefficientlistden)]==0.0 &&
        sum(i->abs(coefficientlistnum[i]-coefficientlistden[i-1]),
        2:length(coefficientlistden))==0.0
    println("All the points are fixed points")
    else
    ffgg=pairoffunctions(coefficientlistnumcomplex, coefficientlistdencomplex)
    Polynumerator=Poly(coefficientlistnumcomplex)
    Polydenominator=Poly(coefficientlistdencomplex)
    Polyx=Poly([0.0+ 0.0* im,1.0+0.0*im])
    lookingzeros= Polynumerator- Polydenominator * Polyx
    fix=roots(lookingzeros)
    le=length(fix)
        fixed=[(complex(fix[i]),1.0+0.0*im) for i in 1:le]
        if abs(ffgg[1](complex(0.0),complex(0.0)))==0.0 &&
            abs(ffgg[2](complex(0.0),complex(0.0)))==0.0
       newfixed=[(complex(0.0),complex(0.0))]
        else
       newfixed=[]
        end
        if abs(ffgg[1](complex(1.0),complex(0.0)))!=0.0 &&
            abs(ffgg[2](complex(1.0),complex(0.0)))==0.0
       morenewfixed=[(complex(1.0),complex(0.0))]
        else
       morenewfixed=[]
        end
        return newfixed, morenewfixed, fixed
    end
end


fixedPointsofaIrreduciblePairex=
    fixedPointsofaIrreduciblePair(coefficientlistnumcompila,
    coefficientlistdencompila,9)[3]

function newstep(hpair::Tuple{Function,Function},
                    iter::Int64,
                    iterprecision::Int64,
                    hpoint::Tuple{Complex{Float64},Complex{Float64}})
    point = hpoint
    tol=1.0/10^(iterprecision)
    number = 0
    imagepoint = rationalFunction(hpair,point,iterprecision)
    while
    (chordalMetric(point, imagepoint,iterprecision) > tol)&& (number < iter)
    point = imagepoint
    imagepoint = rationalFunction(hpair,point,iterprecision)
    number=number+1
    end
    return [imagepoint,number]
end

newstep(hpaircompila, 11, 3, (1.0+0.0*im, 0.0+1.0*im))


function rectangle(xinterval::Tuple{Float64,Float64}=(-1.5,1.5),
                    yinterval::Tuple{Float64,Float64}=(-1.5,1.5),
                    expprecision=10)
    tol=1.0/(2^expprecision)
    a=xinterval[1]
    b=xinterval[2]
    c=yinterval[1]
    d=yinterval[2]
    red=[(complex(r,i),complex(1.0)) for i=d:-tol:c, r=a:tol:b]
    return red
end


function newstep(hpair::Tuple{Function,Function}, iter::Int64,
        iterprecision::Int64,
        hrectangle::Array{Tuple{Complex{Float64},Complex{Float64}},2})
    size1=size(hrectangle)
    result=[newstep(hpair, iter, iterprecision, hrectangle[i,j])
        for i=1:size1[1], j=1:size1[2]]
    return result
end

#Example
rectanglecompila=rectangle((0.0,1.0),(0.0,1.0),1)
newstep(hpaircompila, 11, 3, rectanglecompila)

function positionuptotolerance(fixedPointList::Array{Tuple{Complex{Float64},
                            Complex{Float64}},1},
                            aproxprecision::Int64,
                            twotuple::Tuple{Complex{Float64},Complex{Float64}})
    pos=0
    it=1
    le = length(fixedPointList)
    tol=1/10^(aproxprecision)
    while (it < le+1)
        if (chordalMetric(twotuple, fixedPointList[it],aproxprecision) < tol)
            pos = it
        end
    it=it+1
    end
    return convert(Int64,pos)
end

positionuptotolerance(fixedPointsofaIrreduciblePairex,3,(-0.5 - 0.28867513459481275im,1.0 + 0.0im) )

function position_iteration_upto_tolerances(hpair::Tuple{Function,Function},
        fixedPointList::Array{Tuple{Complex{Float64},Complex{Float64}},1},
        iter::Int64, iterprecision::Int64,
        aproxprecision::Int64,
        twotuple::Tuple{Complex{Float64},Complex{Float64}})
    endpoint_iterations=newstep(hpair, iter, iterprecision, twotuple)
    endpoint=endpoint_iterations[1]
    iterations=convert(Int64,endpoint_iterations[2])
    pos=positionuptotolerance(fixedPointList, aproxprecision, endpoint)
    return (pos,   iterations)
end

position_iteration_upto_tolerances(hpaircompila,fixedPointsofaIrreduciblePairex, 11, 3, 2, (1.0+0.0*im, 0.0+1.0*im))


function position_iteration_upto_tolerances(hpair::Tuple{Function,Function},
        fixedPointList::Array{Tuple{Complex{Float64},Complex{Float64}},1},
        iter::Int64, iterprecision::Int64, aproxprecision::Int64,
        hrectangle::Array{Tuple{Complex{Float64},Complex{Float64}},2})
    size1=size(hrectangle)
    result=[position_iteration_upto_tolerances(hpair,fixedPointList, iter,
        iterprecision, aproxprecision,
        hrectangle[i,j]) for i=1:size1[1], j=1:size1[2]]
    return result
end

position_iteration_upto_tolerances(hpaircompila,fixedPointsofaIrreduciblePairex,11, 3, 2, rectanglecompila)

function fixedPointListex_matrixpositioninterations_RationalFunction(
        coefficientlistnum::Array{T,1},
        coefficientlistden::Array{T,1},
        preexpresolution::Int=8,
        preiteration_max::Int=25,
        preiterprecision::Int64=3,
        preaproxprecision::Int64=3,
        premodeldomain::AbstractString="localrectangle",
        prerectanglesidesdomain::Tuple{Float64,Float64,Float64,Float64}=
        (-1.5,1.5,-1.5,1.5)) where {T<:Number}
    xinterv=(prerectanglesidesdomain[1],prerectanglesidesdomain[2])
    yinterv=(prerectanglesidesdomain[3],prerectanglesidesdomain[4])
    rect=rectangle(xinterv,  yinterv, preexpresolution)
        hpair=pairoffunctions(coefficientlistnum, coefficientlistden)
        coefficientlistnumcomplex=complex(coefficientlistnum)
        coefficientlistdencomplex=complex(coefficientlistden)
    fixedPointListexcomplete=
        fixedPointsofaIrreduciblePair(coefficientlistnumcomplex,
        coefficientlistdencomplex,
        preaproxprecision)
    fixedPointListexcompletea=fixedPointListexcomplete[1]
    if length(fixedPointListexcomplete[2])==0
    fixedPointListexcompleteb=fixedPointListexcompletea
    else
    fixedPointListexcompleteb=union(fixedPointListexcompletea,
        [(complex(1.0),complex(0.0))])
    end
    fixedPointListexample=
    union(fixedPointListexcompleteb,fixedPointListexcomplete[3])
    positer=position_iteration_upto_tolerances(hpair,fixedPointListexample,
    preiteration_max, preiterprecision,
    preaproxprecision, rect)
    return   fixedPointListexample, positer, fixedPointListexcomplete
end

fmf=fixedPointListex_matrixpositioninterations_RationalFunction(
    coefficientlistnumcompila,
    coefficientlistdencompila,3)

function plot_matrixpositer_RationalFunction(
    fixedpointsmatrixpositer,iteration_max::Int=25,
    thecolorstrategy::AbstractString="positionplusiteration",
    themodel::AbstractString="localrectangle",
    therectanglesides::Tuple{Float64,Float64,Float64,Float64}=
    (-1.5,1.5,-1.5,1.5))
    numberoffixedPointListex=length(fixedpointsmatrixpositer[1])
    positer=fixedpointsmatrixpositer[2]
    ab=size(positer)
        if thecolorstrategy=="positionplusiteration"
        numberofcolors=numberoffixedPointListex+1
        integermatrix=
        [(positer[i,j][1]+(1.0-(positer[i,j][2])/iteration_max))+1.0  for i=1:ab[1], j=1:ab[2]];
        integermatrix[1,1]=1;integermatrix[ab[1],ab[2]]=numberofcolors-2
        elseif thecolorstrategy=="iteration"
        numberofcolors=iteration_max+1
        integermatrix=[positer[i,j][2]+0.25 for i=1:size(positer)[1],
            j=1:size(positer)[2]];integermatrix[1,1]=iteration_max
        else thecolorstrategy=="positionfixedpoints"
        numberofcolors=numberoffixedPointListex+1
        integermatrix=[positer[i,j][1]+0.25 for i=1:size(positer)[1],
            j=1:size(positer)[2]]; integermatrix[1,1]=0;integermatrix[ab[1],
            ab[2]]=numberoffixedPointListex
        end
    isinfinity=length(fixedpointsmatrixpositer[3][2])
    if  isinfinity==0
        seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5)]
    else
        seed1=[Colors.RGB(0.0,0.0,0.0),Colors.RGB(0.5,0.5,0.5),Colors.RGB(1.0,1.0,0.0)]
    end
    seed=union(seed1,[Colors.RGB(1.0,0.0,0.0),Colors.RGB(0.0,1.0,0.0),Colors.RGB(0.0,0.0,1.0)])
    gseed=Colors.distinguishable_colors(numberofcolors, seed)
    fpcm=PyPlot.ColorMap(gseed)
    img=PyPlot.imshow(integermatrix, cmap=fpcm, extent=[therectanglesides[1],
        therectanglesides[2],
        therectanglesides[3],therectanglesides[4]])
    return  fixedpointsmatrixpositer[3],gseed, img
end

function aux_plottingBasinsUnivariateRationalFunctions(
    coefficientlistnum::Array{T,1},
    coefficientlistden::Array{T,1},
    expresolution::Int=8,
    iterationmax::Int=25,
    iterprecision::Int64=3,
    aproxprecision::Int64=3;
    colorstrategy::AbstractString="positionplusiteration",
    model::AbstractString="localrectangle",
    rectanglesides::Tuple{Float64,Float64,Float64,Float64}=
        (-1.5,1.5,-1.5,1.5)) where {T<:Number}
    maxdegree=max(length(coefficientlistnum),length(coefficientlistden))
    samedegreecoefficientlistnum=zeros(Complex{Float64}, maxdegree)
    for i in 1:length(coefficientlistnum)
        samedegreecoefficientlistnum[i]=coefficientlistnum[i]
    end
    samedegreecoefficientlistden=zeros(Complex{Float64}, maxdegree)
    for i in 1:length(coefficientlistden)
        samedegreecoefficientlistden[i]=coefficientlistden[i]
    end
    preexpresolution=expresolution
    preiterationmax=iterationmax
    preiterprecision=iterprecision
    preaproxprecision=aproxprecision
    premodeldomain=model
    prerectanglesidesdomain=rectanglesides
    fixedpointsmatrixpairpositerex=
    fixedPointListex_matrixpositioninterations_RationalFunction(
    samedegreecoefficientlistnum,
    samedegreecoefficientlistden, preexpresolution,preiterationmax,
    preiterprecision,preaproxprecision,
     premodeldomain,  prerectanglesidesdomain)
    fixedpointsmatrixpairpositer3seedimg=
    plot_matrixpositer_RationalFunction(fixedpointsmatrixpairpositerex,
                                        iterationmax,
                                        colorstrategy,
                                        model,
                                        rectanglesides)
    return  fixedpointsmatrixpairpositer3seedimg
end

aux_result=aux_plottingBasinsUnivariateRationalFunctions(coefficientlistnumcompila,
    coefficientlistdencompila)

function plottingBasinsUnivariateRationalFunctions(
    coefficientlistnum::Array{T,1},
    coefficientlistden::Array{T,1},
    expresolution::Int=8,
    iterationmax::Int=25,
    iterprecision::Int64=3,
    aproxprecision::Int64=3;
    colorstrategy::AbstractString="positionplusiteration",
    model::AbstractString="localrectangle",
    rectanglesides::Tuple{Float64,Float64,Float64,Float64}=
        (-1.5,1.5,-1.5,1.5)) where {T<:Number}
    maxdegree=max(length(coefficientlistnum),length(coefficientlistden))
    samedegreecoefficientlistnum=zeros(Complex{Float64}, maxdegree)
    for i in 1:length(coefficientlistnum)
        samedegreecoefficientlistnum[i]=coefficientlistnum[i]
    end
    samedegreecoefficientlistden=zeros(Complex{Float64}, maxdegree)
    for i in 1:length(coefficientlistden)
        samedegreecoefficientlistden[i]=coefficientlistden[i]
    end
    preexpresolution=expresolution
    preiterationmax=iterationmax
    preiterprecision=iterprecision
    preaproxprecision=aproxprecision
    premodeldomain=model
    prerectanglesidesdomain=rectanglesides
    fixedpointsmatrixpairpositerex=
    fixedPointListex_matrixpositioninterations_RationalFunction(
    samedegreecoefficientlistnum,
    samedegreecoefficientlistden, preexpresolution,preiterationmax,
    preiterprecision,preaproxprecision,
     premodeldomain,  prerectanglesidesdomain)
    fixedpointsmatrixpairpositer3seedimg=
    plot_matrixpositer_RationalFunction(fixedpointsmatrixpairpositerex,
                                        iterationmax,
                                        colorstrategy,
                                        model,
                                        rectanglesides)
    return  fixedpointsmatrixpairpositer3seedimg[3]
end

result=plottingBasinsUnivariateRationalFunctions(coefficientlistnumcompila,
    coefficientlistdencompila)
end # module
