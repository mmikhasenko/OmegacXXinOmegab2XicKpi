"""
    λ(x,y,z)

Kàllën triangle function
""" 
λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x

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


"""
    pq(m,m1,m2,m3,m0)

Product of the break-up momenta in the three-body decay
    m0 -> m + m3
          ↪ m1 + m2
it gives the phase-space factor x jacobian when projecting to the m-dimension
""" 
pq(m,m1,m2,m3,m0) = m1+m2 < m < m0-m3 ? sqrt(λ(m^2, m1^2, m2^2)*λ(m^2, m0^2, m3^2)) / (2*m*m0) : zero(m)

"""
    Φ2(x,m1,m2)

Two-body phase space factor, x->m1+m2
""" 
Φ2(x,m1,m2) = sqrt(λ(x^2,m1^2,m2^2))/x^2
