function t1(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    return (kₐ+k_d)/(k_d*q) * (1 + k₋/k₊⁺ + k₋*ℓ/p) * (1 - exp(ℓ * p/q))/(1-exp(p/q))
end

function t0(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    if q > p
        p0 = exp(- p*k₋/((q-p)*k₊))
    else 
        p0 = 1
    end
    if q*(k₊)/(k₊+k₋) > p
        return Inf
    else
        t0 = q/(p*k₋) + (1/k₋ + 1/k₊) * (q-p)/p*k₊/k₋+ (1-p0) * ((q-p)/(p/k₊ - (q-p)/k₋)*1/k₋*(1/k₋ + 1/k₊))
        return t0
    end
end

function V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    v_c = (((q/p)^ℓ-1)/((q/p)^(ℓ+1)-1)) * q * k₊⁺/(k₊⁺+k₋)
    v_f = p
    p0 = t1(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1(p,q,kₐ,k₊,k₋,Eₐ,E₊) + t0(p,q,kₐ,k₊,k₋,Eₐ,E₊))
    return v_f * (1-p0) + v_c * p0
end


function C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    return kₐ/(kₐ+k_d) * t1(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1(p,q,kₐ,k₊,k₋,Eₐ,E₊) + t0(p,q,kₐ,k₊,k₋,Eₐ,E₊))
end

function Fₜ₊(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ)
    if q > p
    p0 = exp(- p*k₋/((q-p)*k₊))
    else 
    p0 = 1
    end
    function te(p,q,kₐ,k₊,k₋,Eₐ,E₊)
        if q*(k₊)/(k₊+k₋) > p
            return Inf
        else
            te = (q/(p*k₋) + (1/k₋ + 1/k₊) * (q-p)/p*k₊/k₋) * (exp(-k₋*(ℓₚ-ℓ)/(q-p))) + (1-p0) * ((q-p)/(p/k₊ - (q-p)/k₋)*1/k₋*(1/k₋ + 1/k₊))
            return te
        end
    end
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    k₊⁺ = k₊ * exp(E₊)
    return 1 -  te(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1(p,q,kₐ,k₊,k₋,Eₐ,E₊)+t0(p,q,kₐ,k₊,k₋,Eₐ,E₊))
end
function C₀(Eₐ)
    local kₐ = 1
    k_d = kₐ * exp(-Eₐ)
    return kₐ/(kₐ+k_d)
end