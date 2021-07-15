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
    legendfontsize=9, legend =:topright)
#
using CSV
using JSON
using Parameters
using DataFrames
using DelimitedFiles

@userplot XicKpiFitPlot
@recipe function f(hp::XicKpiFitPlot; binning=range(extrema(hp.args[2])...,length=84))
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
    labels = [
        L"\mathit{\Omega}_b^-\rightarrow\mathit{\Xi}_c^+ K^-\pi^-",
        L"\mathrm{Comb. background}",
        L"\mathit{\Omega}_b^-\rightarrow\mathit{\Xi}_c^{+}K^-\rho^-(\rightarrow\pi^-\pi^0)",
        L"\mathit{\Omega}_b^-\rightarrow\mathit{\Xi}_c^{+\prime}(\rightarrow\mathit{\Xi}_c^{+}\gamma\,)\,K^-\pi^-",
        L"\mathit{\Omega}_b^-\rightarrow\mathit{\Xi}_c^+ K^-K^-"]
    for i in 1:5
        @series begin
            label := labels[i]
            seriescolor := i
            (model[i], w, 300)
        end
    end
    xlims --> extrema(binning)
    ylims --> (0,:auto)
    # 
    bin_width = Int(round(1e3*(binning[2]-binning[1])))
    xguide := L_ΞcKπ
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

@unpack feeddown = settings["XicKpi_fit"]
XicKpi_fitrange = Tuple(settings["XicKpi_fit"]["XicKpi_fitrange"])
@unpack f_narrow_gauss, σ2_over_σ1, alpha_bgd, fKK2Kpi = settings["XicKpi_fit"]

# read data
# data = CSV.File(joinpath("data", "tnominal.csv")) |> DataFrame |> d->select(d, :mΞcKπ)
# CSV.write(joinpath("data", "XicKpi_fit.csv"), data)
# 
data = CSV.File(joinpath("data", "XicKpi_fit.csv")) |> DataFrame
Nev = size(data,1)

signal = Normalized(
    FDoubleGaussFixedRatio(
        Ext(mΩb=6.045, σ=0.013, r=f_narrow_gauss, n=σ2_over_σ1)),
    XicKpi_fitrange) |> x->fixpars(x,(:r,:n))
# 
plot(signal)
# 
combbgd = Normalized(FExp(Ext(τ=-5.1,)), XicKpi_fitrange) |> x->fixpar(x,:τ)
plot(combbgd)

# feeddown
function pdffromXY(XY, range)
    f = FTabulated(XY[:,1], XY[:,2])
    return Normalized(f, range)
end
tXi = pdffromXY(readdlm(feeddown["Xicp"])[:,[1,2]], XicKpi_fitrange)
tRh = pdffromXY(readdlm(feeddown["rho"] )[:,[1,2]], XicKpi_fitrange)
tKK = pdffromXY(readdlm(feeddown["KK"]  )[:,[1,2]], XicKpi_fitrange)
# 
model1 = FSum(
    [signal, combbgd, tRh, tXi, tKK],
    Ext(fS=0.6Nev, fB=0.2Nev, cfRh=0.1Nev, cfXi=0.1Nev, cfKK=0.02Nev))
# 
xickpifitplot(model1, data.mΞcKπ)

# do two tries:
# first
fit1 = fit(Extended(NegativeLogLikelihood(model1, data.mΞcKπ)))
# 
xickpifitplot(fit1.best_model, data.mΞcKπ)

# fix fraction of fKK
fKK = fit1.parameters.fS * fKK2Kpi
model2 = updatepars(model1, fit1.parameters) |> x->fixpar(model1, :cfKK, fKK)

# second
fit2 = fit(Extended(NegativeLogLikelihood(model2, data.mΞcKπ)))
# 
xickpifitplot(fit2.best_model, data.mΞcKπ)

savefig(joinpath("plots", "default", "XicKpi_fit.pdf"))

# determine the cut range
twoσrange = let 
    μ = fit2.parameters.mΩb
    σ1 = fit2.parameters.σ
    σ2 = σ1*σ2_over_σ1
    f = f_narrow_gauss
    # 
    σ = sqrt(f*σ1^2+(1-f)*σ2^2)
    μ .+ (2σ) .* (-1,1)
end
# 

# determine the yield from the integral
bgd_integral = let
    bgd_fraction = integral(fit2.best_model[2], twoσrange)
    bgd_fraction *
        fit2.measurements.fB / fit2.parameters.fB
end

# full yield
signal_integral = let
    signal_fraction = integral(fit2.best_model[1], twoσrange)
    signal_fraction *
        fit2.measurements.fS / fit2.parameters.fS
end
signal_selected = size(filter(x->inrange(x.mΞcKπ,twoσrange), data), 1)


# select signal
writejson(joinpath("results", "default", "XicKpi_fit.json"),
    transformdictrecursively!(Dict(
        :fit => Dict(
            :fixedpars => fixedpars(model2),
            :freepars => fit2.measurements),
        :inferred => Dict(
            :signal_mass => fit2.measurements.mΩb*1e3,
            :signal_integral => signal_integral,
            :signal_selected => signal_selected,
            :twoσrange => round.(twoσrange, digits=4),
            :bgd_integral => bgd_integral)
    ), ifmeasurementgivestring)
)
#
