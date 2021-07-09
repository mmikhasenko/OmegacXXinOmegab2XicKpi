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
using Interpolations
using Statistics
using DataFrames
using DelimitedFiles

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
data = CSV.File(joinpath("data", "tnominal.csv")) |> DataFrame

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
model0 = MixedModel(
    [combbgd, signal, tRh, tKK, tXi],
    Ext(fB=0.2, fS=0.6, cfRh=0.1, cfKK=0.02))
plot(data.mΞcKπ, model0)

# do two tries:
# first
fit1 = fit(model1, data.mΞcKπ)

# fix fraction of fKK
fKK = fractionvalues(fit1.best_model)[2] * fKK2Kpi
model2 = updatepars(model1, fit1.parameters) |> x->fixpar(model1, :cfKK, fKK)

# second
fit2 = fit(model2, data.mΞcKπ)

let bins = range(XicKpi_fitrange..., length=44)
    best_model = fit2.best_model
    Ns = scaletobinneddata(length(tnom.mΞcKπ), bins)
    plot(best_model.components[4], fractionvalues(best_model)[4]*Ns, lw=1.5, lab=L"\mathit{\Omega}_b^-\rightarrow\mathit{\Xi}_c^+ K^-K^-",c=6)
    plot!(xlab=L_ΞcKπ, ylab=L_cand_ofXeV(Int(round(1e3*(bins[2]-bins[1])))), lw=1.5, c=2)
    plot!(best_model.components[1], fractionvalues(best_model)[1]*Ns, lw=1.5, lab=L"\mathrm{Comb. background}")
    plot!(best_model.components[5], fractionvalues(best_model)[5]*Ns, lw=1.5, lab=L"\mathit{\Omega}_b^-\rightarrow\mathit{\Xi}_c^{+\prime}(\rightarrow\mathit{\Xi}_c^{+}\gamma)\,K^-\pi^-")
    plot!(best_model.components[3], fractionvalues(best_model)[3]*Ns, lw=1.5, lab=L"\mathit{\Omega}_b^-\rightarrow\mathit{\Xi}_c^{+}K^-\rho^-(\rightarrow\pi^-\pi^0)")
    plot!(tnom.mΞcKπ, best_model, lab=L"\mathrm{Model}", datalabel=L"\mathrm{Data}", c=1, lc=:red, lw=1)
    plot!(xlims=XicKpi_fitrange, ylims=(0,:auto))
end
savefig(joinpath("plots", "default", "XicKpi_fit.pdf"))

# determine the cut range
twoσrange = let 
    μ = fit_nominal.parameters.mΩb
    σ1 = fit_nominal.parameters.σ
    σ2 = σ1*σ2_over_σ1
    f = f_narrow_gauss
    # 
    σ = sqrt(f*σ1^2+(1-f)*σ2^2)
    μ .+ (2σ) .* (-1,1)
end
# 

# determine the yield from the integral
bgd_integral = let
    bgd_fraction = integrals(fit_nominal.best_model, twoσrange)[1]
    bgd_fraction * length(data.mΞcKπ) *
        fit_nominal.measurements.fB / fit_nominal.parameters.fB
end

# full yield
signal_integral = let
    signal_fraction = integrals(fit_nominal.best_model, twoσrange)[2]
    signal_fraction * length(data.mΞcKπ) *
        fit_nominal.measurements.fS / fit_nominal.parameters.fS
end
signal_selected = size(filter(x->inrange(x.mΞcKπ,twoσrange), data), 1)


# select signal
writejson(joinpath("results", "default", "XicKpi_fit.json"),
    transformdictrecursively!(Dict(
        :fit => Dict(
            :fixedpars => fixedpars(model),
            :freepars => fit_nominal.measurements),
        :inferred => Dict(
            :signal_mass => fit_nominal.measurements.mΩb*1e3,
            :signal_integral => signal_integral,
            :signal_selected => signal_selected,
            :twoσrange => round.(twoσrange, digits=4),
            :bgd_integral => bgd_integral)
    ), ifmeasurementgivestring)
)
#
