
@with_kw struct stateHc{T}
    n::Int
    m::T
    Γ::T
end

const mΩc3000 = 3000.4e-3; const mΩc3000sq = mΩc3000^2
const mΩc3050 = 3050.2e-3; const mΩc3050sq = mΩc3050^2
const mΩc3065 = 3065.6e-3; const mΩc3065sq = mΩc3065^2
const mΩc3090 = 3090.2e-3; const mΩc3090sq = mΩc3090^2
const mΩc3119 = 3119.1e-3; const mΩc3119sq = mΩc3119^2
const mΩc3188 = 3188.0e-3; const mΩc3188sq = mΩc3188^2
#
const ΓΩc3000 =  4.5e-3  # 2998.5 ± 0.9   7.0 ± 2.0   44.5 ± 7.5
const ΓΩc3050 =  0.8e-3  # 3050.2 ± 0.3   0.7 ± 0.5   40.1 ± 6.8
const ΓΩc3065 =  3.5e-3  # 3065.9 ± 0.4   2.3 ± 0.8   51.8 ± 8.2
const ΓΩc3090 =  8.7e-3  # 3091.1 ± 1.1   9.4 ± 3.0   49.3 ± 10.1
const ΓΩc3119 =  1.1e-3  # 3091.1 ± 1.1   9.4 ± 3.0   49.3 ± 10.1
const ΓΩc3188 = 60.0e-3  #
#
statesΩc = [let
    ms = Symbol("mΩc$n"); m = @eval $ms
    Γs = Symbol("ΓΩc$n"); Γ = @eval $Γs
    stateHc(n=n, m=m, Γ=Γ)
end for n in [3000, 3050, 3065, 3090, 3119, 3188]]
#
statesΩc4 = statesΩc[1:4]
#
