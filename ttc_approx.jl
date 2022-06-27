function V̄(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    k_d = kₐ * exp(-Eₐ)
    k₊⁺ = k₊ * exp(E₊)
    return q*((q/p)^ℓ-1)/((q/p)^(ℓ+1)-1) * (k₊⁺/(k₊⁺+ k₋))
end
using SpecialFunctions
function EFₜ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ)
    k₊⁺ = k₊ * exp(E₊)
    return 1 - 1/2 * (1 - erf((ℓₚ-ℓ-(q-p)/k₋)/√(2*(p+q)/k₋))) / (1+ p*k₋*(k₊⁺+k₋)/(q*k₊⁺)*(1+exp(Eₐ)))
end