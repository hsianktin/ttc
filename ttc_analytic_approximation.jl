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
# function C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊)
#     k₊⁺ = k₊ * exp(E₊)
#     k_d = kₐ * exp(-Eₐ)
#     t1 = (kₐ+k_d)/(k_d*(q-p)) * (1 + k₋/k₊⁺)
#     t0 = (q/(k₋ * p) + 1/kₐ)/  (1-k₊/(k₊ + p/(q-p) *k₋))
#     return kₐ/(kₐ+k_d) * t1/(t1 + t0)
# end


## Testing about Coupling coefficient
q = 30
kₐ = 100
Eₐ = 3
k₊ = 0.3
k₋ = 0.4
L = 335
E₊ = 2
ℓ = 4
using Plots
plotly()
p = [i for i in 2:0.2:40]
Cₑₛₜ = [C₊(i,q,kₐ,k₊,k₋,Eₐ,E₊) for i in p]

plot(p,Cₑₛₜ)

CSV.write("fig/plot_df_C_p.csv", DataFrame(p=p,eta=p./(q*k₊/(k₊+k₋)),C=Cₑₛₜ))

using DataFrames
using CSV
df = CSV.read("fig/simu_df_ps_prediction.csv",DataFrame)
plot!(df.p, df.C, label="L=335", line=:scatter)
xlabel!("p")
ylabel!("C")

df = CSV.read("fig/simu_df_ps_long_L.csv",DataFrame)
plot!(df.p, df.C, label="L→∞", line=:scatter)
xlabel!("p")
ylabel!("C")

plot!(p, [C₀(Eₐ) for i in p], label="C₀", line=:dash)

Eₐs = [i for i in -2:0.25:7];
p = 20;
q = 30;
E₊ = 2;
Cₑₛₜ = [C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊) for Eₐ in Eₐs]
plot(Eₐs,Cₑₛₜ)
df = CSV.read("fig/simu_df_Ebs_2.csv",DataFrame)
plot!(df.Eᵦ, df.C, label="L=335", line=:scatter)
xlabel!("Eₐ")
ylabel!("C")

E₊s = [i for i in 0:0.1:5]
p = 20;
q = 30;
Eᵦ = 3;
Cₑₛₜ = [C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊) for E₊ in E₊s]
plot(E₊s,Cₑₛₜ)
df = CSV.read("fig/simu_df_E_cs.csv",DataFrame)
plot!(df.E_c, df.C, label="L=335", line=:scatter)

## Testing about protection fraction
q = 30
kₐ = 100
Eₐ = 3
k₊ = 0.3
k₋ = 0.4
L = 335
E₊ = 2
ℓ = 4
ℓₚ = 27
ps = [i for i in 2:0.2:30]
Fₑₛₜ =  [Fₜ₊(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ) for p in ps] |> # replace nan with 0
    x -> [isnan(i) ? 0 : i for i in x]

plot(ps,Fₑₛₜ)
CSV.write("fig/plot_df_F_p.csv", DataFrame(p=ps,eta=ps./(q*k₊/(k₊+k₋)),F_T=Fₑₛₜ))

using DataFrames
using CSV
df = CSV.read("fig/simu_df_ps_prediction.csv",DataFrame)
plot!(df.p, df.fractions_T_protected, label="L=335", line=:scatter)
xlabel!("p")
ylabel!("Fₜ")

df = CSV.read("fig/simu_df_ps_long_L.csv",DataFrame)
plot!(df.p, df.fractions_T_protected, label="L→∞", line=:scatter)


Eₐs = [i for i in -2:0.25:7];
p = 20;
q = 30;
E₊ = 2;
Fₑₛₜ = [Fₜ₊(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ) for Eₐ in Eₐs]
plot(Eₐs,Fₑₛₜ)
df = CSV.read("fig/simu_df_Ebs_2.csv",DataFrame)
plot!(df.Eᵦ, df.fractions_T_protected, label="L=335", line=:scatter)
xlabel!("Eₐ")
ylabel!("F_T")
df = CSV.read("fig/simu_df_Ebs_L_inf_2.csv",DataFrame)
plot!(df.Eᵦ, df.fractions_T_protected, label="L→∞", line=:scatter)

E₊s = [i for i in 0:0.1:5]
p = 20;
q = 30;
Eᵦ = 3;
Fₑₛₜ = [Fₜ₊(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ) for E₊ in E₊s]
plot(E₊s,Fₑₛₜ)
df = CSV.read("fig/simu_df_E_cs.csv",DataFrame)
plot!(df.E_c, df.fractions_T_protected, label="L=335", line=:scatter)
df = CSV.read("fig/simu_df_E_cs_L_inf.csv",DataFrame)

xlabel!("E₊")
ylabel!("F_T")


## Testing about the coupled velocity
q = 30
kₐ = 100
Eₐ = 3
k₊ = 0.3
k₋ = 0.4
L = 335
E₊ = 2
ℓ = 4
using Plots
plotly()
ps = [i for i in 2:0.2:29]
Vₑₛₜ = [V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for p in ps]
plot(ps,Vₑₛₜ)

using DataFrames
using CSV
df = CSV.read("fig/simu_df_ps_prediction.csv",DataFrame)
plot!(df.p, df.v_eff_translation, label="L=335", line=:scatter)
xlabel!("p")
ylabel!("V_translation")

df = CSV.read("fig/simu_df_ps_long_L.csv",DataFrame)
plot!(df.p, df.v_eff_translation, label="L→∞", line=:scatter)


Eₐs = [i for i in -2:0.25:7];
p = 20;
q = 30;
E₊ = 2;
Vₑₛₜ = [V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for Eₐ in Eₐs]
plot(Eₐs,Vₑₛₜ)
df = CSV.read("fig/simu_df_Ebs_2.csv",DataFrame)
plot!(df.Eᵦ, df.v_eff_translation, label="L=335", line=:scatter)
xlabel!("Eₐ")
ylabel!("C")

E₊s = [i for i in 0:0.1:5]
p = 20;
q = 30;
Eᵦ = 3;
Vₑₛₜ = [V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for E₊ in E₊s]
plot(E₊s,Vₑₛₜ)
df = CSV.read("fig/simu_df_E_cs.csv",DataFrame)
plot!(df.E_c, df.v_eff_translation, label="L=335", line=:scatter)

CSV.write("fig/plot_df_V_E_cs.csv", DataFrame(Ec=E₊s, V=Vₑₛₜ))

E₊s = [i for i in 0:0.5:2]
ps = [i for i in 13:0.2:60]
Vₑₛₜ = [V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for p in ps, E₊ in E₊s]
plot(ps,Vₑₛₜ)
df = CSV.read("fig/merged_df_ps_var_E_cs.csv",DataFrame)
for E₊ in E₊s
    plot!(df.p[df.Ec .== E₊], df.v_eff_transcriptions[df.Ec .== E₊], label="E₊=$E₊", line=:scatter)
end
xlabel!("p")

plot_df = DataFrame(
    p = Float64[],
    V = Float64[],
    Ec = Float64[]
)
for E₊ in E₊s
    for p in ps
        push!(plot_df, [p, V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ), E₊])
    end
end
CSV.write("fig/plot_df_V_p_E_cs.csv", plot_df)