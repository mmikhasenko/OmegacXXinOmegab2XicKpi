using Plots: length
cd(joinpath(@__DIR__,".."))
using Pkg
Pkg.activate(".")
# 
using OmegacXXinOmegab2XicKpi
# 
using AlgebraPDF

using Plots
using Plots.PlotMeasures
using LaTeXStrings
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))
#
using CSV
using JSON
using Parameters
using DataFrames

# minuit stuff: fixing parameters
"""
    limit_argument_key(parameters, symbol_range)

Key argument for parameter limit in minuit.
```jldoctest
julia> limit_argument_key((a=1.1, b=2.2, c=3.3), a => (1,2))
(limit_x0 = (1,2), )
```
"""
function limit_argument_key(parameters, symbol_range)
    symbol, range = symbol_range
    key = Symbol("limit_x", find_index_xi(parameters, symbol)-1)
    NamedTuple{(key,)}((range,))
end
find_index_xi(parameters, symbol) = findfirst(x->x==symbol, collect(keys(parameters)))

# 
readresults(filename) = transformdictrecursively!(
    readjson(joinpath("results", "default", filename)),
    ifstringgivemeasurement)
# 
e2m(e) = e/1e3+(mΞc+mK)
m2e(m) = (m-(mΞc+mK))*1e3

phsp_e(e) = pq(e2m(e),mΞc,mK,mπ,mΩb)

#            _|              _|      
#  _|_|_|    _|    _|_|    _|_|_|_|  
#  _|    _|  _|  _|    _|    _|      
#  _|    _|  _|  _|    _|    _|      
#  _|_|_|    _|    _|_|        _|_|  
#  _|                                
#  _|                                

@userplot XicKFitPlot
@recipe function f(hp::XicKFitPlot; binning=range(extrema(hp.args[2])...,length=84))
    model, xdata = hp.args
    @series begin
        seriestype := :scatter
        label := L"\mathrm{Data}"
        markershape := :+
        ([-2.0], [-1.0])
    end
    # lines
    w = scaletobinneddata(binning)
    @series begin
        label := L"\mathrm{Model}"
        linewidth --> 2
        linecolor --> :red
        (model, w, 300)
    end
    Nterms = length(model)
    @series begin
        label := L"\mathrm{Combinatirial\,\,bgd}"
        seriescolor := 3
        (model[Nterms], w, 300)
    end
    @series begin
        label := L"\mathrm{Non}\textrm{-}\mathrm{resonant\,\,bgd.}"
        seriescolor := 2
        (model[Nterms-1], w, 300)
    end
    # (+model[Nterms-1]
    xlims --> extrema(binning)
    ylims --> (0,:auto)
    # 
    bin_width = round(binning[2]-binning[1]; digits=1)
    xguide := L_ΞcK_Ξc_K
    yguide := L_cand_ofXeV(bin_width, "MeV")
    #
    seriestype := :stephist
    bins := binning
    label := "" 
    seriescolor := :black
    (xdata,)
end


#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  


settings = readjson("settings.json")
# 
XicK_fitrange = Tuple(settings["XicK_fit"]["XicK_fitrange"])

# select around Omegab peak
twoσrange = Tuple(readresults("XicKpi_fit.json")["inferred"]["twoσrange"])

# read data
data = CSV.File(joinpath("data", "tnominal.csv")) |> DataFrame |> 
    d->filter(x->inrange(x.mΞcKπ, twoσrange), d)
Nev = size(data,1)

# get parameters
@unpack blattweisskopfR = settings["XicK_fit"]
@unpack orbital_momentum = settings["XicK_fit"]

bgd_freepars = d2nt(
        readresults("background_fit_results.json")["massfit"]["freepars"];
    process=x->x.val)

@unpack bgd_integral = readresults("XicKpi_fit.json")["inferred"]

resolutions, resolutions_parameters = let
    b = readresults("resolution.json")
    b["massfit"]["OmegacXXsigma"]..1, 
        b["massfit"]["parameters"]
end
resolutions_parameters

resolution_model =
    FunctionWithParameters(
        (Δm;p)->(Δm>0 ? Δm/1e3*(p.p1+p.p2*Δm/1e3) : 0.0)+1e-6; # function itself
        p=Ext(p1=resolutions_parameters["p1"][1],
              p2=resolutions_parameters["p2"][1])) |>
    d->fixpars(d,(:p1,:p2))
#
plot(resolution_model, 0, 220, lab="",
    xlab=L_ΞcK_Ξc_K, ylab=L"\mathrm{resolution\,\,[MeV]}")

#########################################################################
# components

struct OmegacBW{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
function func(bw::OmegacBW, x::NumberOrTuple; p=pars(bw))
    Δm,Γ = (getproperty(p,s) for s in keys(bw.p))
    amplitudeBWenergydep(e2m(x), e2m(Δm), Γ/1e3; m1=mΞc, m2=mK, L=0, R=1.5)
end

plot(abs2(OmegacBW((Δm1=100,Γ1=20.3))), XicK_fitrange...)


phsp = FunctionWithParameters((x;p)->phsp_e(x), ∅)
plot(phsp, XicK_fitrange..., lab="", xlab=L_ΞcK_Ξc_K)

signalpdfs = []
for i in 1:5
    #
    varnames = (Symbol("Δm"*string(i)), Symbol("Γ"*string(i)))
    Ωc = statesΩc[i] # default values
    nt_mΓ = NamedTuple{varnames}((m2e(Ωc.m), Ωc.Γ*1e3))
    #
    sigi = Normalized(abs2(OmegacBW(nt_mΓ))*phsp, XicK_fitrange)
    # 
    sigixG = convGauss(sigi, resolution_model(m2e(Ωc.m))) # resolution_model resolutions[i]
    push!(signalpdfs, sigixG)
end
@time plot(); plot!.(signalpdfs); plot!()

#  - BW S-wave
@makefuntype OmegacBWg(x;p) = 1/(e2m(p.Δm0)^2 - e2m(x)^2 - 1im*p.g^2*Φ2(e2m(x)+1e-6im,mΞc,mK))
#
signalpdf0 = Normalized(
        abs2(OmegacBWg((Δm0=4.0, g=0.5)))*phsp,
        XicK_fitrange)# |>
    # d->convGauss(d, resolution_model)
plot(signalpdf0)

# backgrounds
#  - nonresonant
nonresbgd = Normalized(phsp, XicK_fitrange)
plot(nonresbgd)

#  - combinatprial background
combbgd = Normalized(FPowExp(Ext(α=0.4, β=-5e-3)), XicK_fitrange) |>
    d->fixpars(d, bgd_freepars)
# 
plot(combbgd)

#########################################################################
# full model

model_nominal = 
    FSum([signalpdf0, signalpdfs[1:4]..., combbgd, nonresbgd],
        Ext(f0=0.05Nev, f1=0.15Nev, f2=0.15Nev, f3=0.25Nev, f4=0.15Nev, fb1=0.02Nev, fb2=0.02Nev)) |>
    d->fixpar(d, :fb1, bgd_integral.val)
# 
@time xickfitplot(model_nominal, data.eΞcK)
#########################################################################
# fit
@time fit_nominal = fit(
    Extended(NegativeLogLikelihood(model_nominal, data.eΞcK));
    limit_argument_key(freepars(model_nominal), :Γ2 => (1e-3,3.0))...)
#
xickfitplot(fit_nominal.best_model, data.eΞcK)
savefig(joinpath("plots", "default", "XicK_fit_unconstrained.pdf"))

#########################################################################
# save

writejson(joinpath("results", "default", "XicK_fit.json"), transformdictrecursively!(
    Dict(
        :unconstrained_fit => Dict(
            :fixedpars => fixedpars(model_nominal),
            :freepars => totuple(fit_nominal.measurements)
            ),
    ),
    ifmeasurementgivestring)
)
