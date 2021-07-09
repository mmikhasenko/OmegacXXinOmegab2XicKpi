"""
    BlattWeisskopfFF(pRsq, L::Int)

Commonly used form-factor to regularize high-energy growth of the amplitude in `L`-wave, with `pRsq = (pR)²`
""" 
function BlattWeisskopfFF(z,L::Int)
    (L > 2 || L < 0) && error("Only 0 ≤ L ≤ 2 are codded! L=$L is called.")
    L==0 && return one(z)
    L==1 && return z/(1+z)
    return z^2/(9+3z+z^2)
end

"""
    amplitudeBWenergydep(x, m, Γ; m1, m2, L::Int=0, R=1.5)

Relativistic Breit-Wigner amplitude describing an isolated elastic resonance is a system of two particles with masses `m₁` and `m₂`
in L-wave. A Blatt-Weisskopf form-factor is used to regularize the high-energy behavior. 
""" 
function amplitudeBWenergydep(x, m, Γ; m1, m2, L::Int=0, R=1.5)
    λK = λ(x^2,m1^2,m2^2)
    λK < 0 && return zero(x)+0im
    # 
    p, p0 = sqrt(λK)/(2x), sqrt(λ(m^2,m1^2,m2^2))/(2m)
    z, z0 = (p*R)^2, (p0*R)^2
    FF, FF0 = BlattWeisskopfFF(z,L), BlattWeisskopfFF(z0,L)
    # 
    Γdep = Γ*p/p0*m/x*FF/FF0
    return m*Γ*sqrt(FF/FF0)/(m^2-x^2-1im*m*Γdep)
end
