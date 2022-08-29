function V̄(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    k₊⁺ = k₊ * exp(E₊)
    return q*((q/p)^ℓ-1)/((q/p)^(ℓ+1)-1) * (k₊⁺/(k₊⁺+ k₋))
end

function V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    f1 = (kₐ+k_d)/(k_d*(q-p)) * (1 + k₋*(1/k₊⁺ + ℓ/(p) + ℓ/(p-q)))
    f2 = q/(k₋ * p) + 1/kₐ
    return q*((((q/p)^ℓ-1)/((q/p)^(ℓ+1)-1))* f1/(1 + k₋*(1/k₊⁺ + ℓ/(p) + ℓ/(p-q))) + 1/k₋)/(f1 + f2)
end

using SpecialFunctions
function 𝔼Fₜ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    f1 = (kₐ+k_d)/(k_d*(q-p)) * (1 + k₋*(1/k₊⁺ + ℓ/(p) + ℓ/(p-q)))
    f2 = q/(k₋ * p) + 1/kₐ
    k₊⁺ = k₊ * exp(E₊)
    return 1 - (maximum([((q-p)/k₋-ℓ)/((q-p)/k₋),0]))*1/2 * (1 - erf((ℓₚ-ℓ-(q-p)/k₋)/√(2*(p+q)/k₋))) * f2 / (f1 + f2)
end

function C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    f1 = (kₐ+k_d)/(k_d*(q-p)) * (1 + k₋/k₊⁺)
    f2 = q/(k₋ * p) + 1/kₐ
    return kₐ/(kₐ+k_d) * f1/(f1 + f2)
end

function Prₙ(x)
    # distribution function of standard normal distribution
    return (1 + erf(x/√2)) / 2
end

function Cₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,α,L)
    q̄ = k₊/(k₊ + k₋) * q
    if L/q̄ - 1 /α < 0
        return 0
    else
        Δ = (L*(p/q̄ - 1) - p/α)/√((L/q̄ - 1/α))
        return Prₙ(Δ) * C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    end
end

