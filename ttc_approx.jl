function V̄(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    k₊⁺ = k₊ * exp(E₊)
    return q*((q/p)^ℓ-1)/((q/p)^(ℓ+1)-1) * (k₊⁺/(k₊⁺+ k₋))
end

function V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    t1 = (kₐ+k_d)/(k_d*(q-p)) * (1 + k₋*(1/k₊⁺ + ℓ/(p) + ℓ/(p-q)))
    t2 = q/(k₋ * p) + 1/kₐ
    return q*((((q/p)^ℓ-1)/((q/p)^(ℓ+1)-1))* t1/(1 + k₋*(1/k₊⁺ + ℓ/(p) + ℓ/(p-q))) + 1/k₋)/(t1 + t2)
end

using SpecialFunctions
function 𝔼Fₜ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    t1 = (kₐ+k_d)/(k_d*(q-p)) * (1 + k₋*(1/k₊⁺ + ℓ/(p) + ℓ/(p-q)))
    t2 = q/(k₋ * p) + 1/kₐ
    k₊⁺ = k₊ * exp(E₊)
    return 1 - (maximum([((q-p)/k₋-ℓ)/((q-p)/k₋),0]))*1/2 * (1 - erf((ℓₚ-ℓ-(q-p)/k₋)/√(2*(p+q)/k₋))) * t2 / (t1 + t2)
end
function Fₜ₊(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    t1 = (kₐ+k_d)/(k_d*(q-p)) * (1 + k₋*(1/k₊⁺ + ℓ/(p) + ℓ/(p-q)))
    t2 = q/(k₋ * p) + 1/kₐ
    k₊⁺ = k₊ * exp(E₊)
    return 1 - t2*exp(-k₋*(ℓₚ-ℓ)/(q-p))/(t1+t2)
end


function C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    t1 = (kₐ+k_d)/(k_d*(q-p)) * (1 + k₋/k₊⁺)
    t2 = q/(k₋ * p) + 1/kₐ
    return kₐ/(kₐ+k_d) * t1/(t1 + t2)
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
        Δ = (L*(p/q̄ - 1) - q̄/α)/√((L/q̄ - 1/α))
        return Prₙ(Δ) * C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    end
end

