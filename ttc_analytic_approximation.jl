function C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    t1 = (kₐ+k_d)/(k_d*(q-p)) * (1 + k₋/k₊⁺)
    t2 = (q/(k₋ * p) + 1/kₐ)/  (k₋ /(q/p*k₊ + k₋))
    return kₐ/(kₐ+k_d) * t1/(t1 + t2)
end

q = 30
kₐ = 100
Eₐ = 3
k₊ = 0.4
k₋ = 0.3
L = 335
E₊ = 2
using Plots
plotly()
p = [i for i in 2:0.2:25]
Cₑₛₜ = [C₊(i,q,kₐ,k₊,k₋,Eₐ,E₊) for i in p]
plot(p,Cₑₛₜ)