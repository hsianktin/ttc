function t1(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    # return (kₐ+k_d)/(k_d*q) * (1 + k₋/k₊⁺ + k₋*ℓ/p)  * (1 - exp(ℓ * p/q))/(1-exp(p/q))
    return t1_em(p,q,kₐ,k₊,k₋,Eₐ,E₊)
end
function t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    # return (kₐ+k_d)/(k_d*q) * (1 + k₋/k₊⁺ + k₋*ℓ/p)  * (1 - exp(ℓ * p/q))/(1-exp(p/q))
    return (kₐ+k_d)/(k_d) * (1 + k₋/k₊⁺ + k₋*ℓ/p)  * ( Tₙ(p,q)  ) *((q+p)/q)
end
function t1_em(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    return (kₐ+k_d)/(k_d*q) * (1 + k₋/k₊⁺ + k₋*ℓ/p)  * (1 - exp(ℓ *(p/q)))/(1-exp(p/q))
    # return (kₐ+k_d)/(k_d) * (1 + k₋/k₊⁺ + k₋*ℓ/p)  * 1/(-λₙ(p,q))
end
using LinearAlgebra

function λₙ(p,q)
    M = zeros(ℓ,ℓ)
    for i in 1:ℓ
        M[i,i] = - (p+q)
        if i < ℓ
            M[i,i+1] = q
        end
        if i > 1
            M[i,i-1] = p
        end
    end
    M[ℓ,ℓ] = -q
    return eigmax(M)   
end

# function Tₙ(ℓ)
#     @vars p,q
#     M = zeros(Sym,ℓ,ℓ)
#     for i in 1:ℓ
#         M[i,i] = 1
#         if i < ℓ
#             M[i,i+1] = -q/(p+q)
#         end
#         if i > 1
#             M[i,i-1] = -p/(p+q)
#         end
#     end
#     M[1,1] = 1
#     M[1,2] = -1
#     v = [1/q ; [1/(p+q) for i in 1:ℓ-1]]
#     return simplify((M^(-1)*v)[1])
# end

function Tₙ(p,q,ℓ)
    M = zeros(ℓ,ℓ)
    for i in 1:ℓ
        M[i,i] = 1
        if i < ℓ
            M[i,i+1] = -q/(p+q)
        end
        if i > 1
            M[i,i-1] = -p/(p+q)
        end
    end
    M[1,1] = 1
    M[1,2] = -1
    v = [1/q ; [1/(p+q) for i in 1:ℓ-1]]
    return ((M^(-1)*v)[1])
end

function Tₙ(p,q)
    factor =0
    for k in 1:ℓ
        factor += k*(q/p)^(k-1)
    end
    return   p^(ℓ-1)/q^(ℓ) * factor
end

function t0(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    return t0_em(p,q,kₐ,k₊,k₋,Eₐ,E₊)
end


function t0_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    if q > p
        p0 = (p/k₊) / ((q-p)/k₋ + p/k₊)
    else 
        p0 = 1
    end
    if q*(k₊)/(k₊+k₋) > p
        return Inf
    elseif p < q
        t0 = (q/(p*k₋) + (1/k₋ + 1/k₊) * maximum([(q-p),0])/p*k₊/k₋+ (1-p0)*((q-p)/(p/k₊ - (q-p)/k₋)*1/k₋*(1/k₋ + 1/k₊)))
        return t0
    else
        t0 = 1/k₋
    end
end

function t0_em(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    if q > p
        p0 = 1-exp(- p*k₋/((q-p)*k₊))
    else 
        p0 = 1
    end
    if q*(k₊)/(k₊+k₋) > p
        return Inf
    else
        t0 = (q/(p*k₋) +(1/k₋ + 1/k₊) * (q-p)/p*k₊/k₋+ (p0) * ((q-p)/(p/k₊ - (q-p)/k₋)*1/k₋*(1/k₋ + 1/k₊)))
        return t0
    end
end

function V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    if p != q
    v_c = (((q/p)^ℓ-1)/((q/p)^(ℓ+1)-1)) * q * k₊⁺/(k₊⁺+k₋)
    else
        v_c = (ℓ)/(ℓ+1) * q * k₊⁺/(k₊⁺+k₋)
    end
    v_f = p
    p2 = t1(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1(p,q,kₐ,k₊,k₋,Eₐ,E₊) + t0(p,q,kₐ,k₊,k₋,Eₐ,E₊))
    return v_f * (1-p2) + v_c * p2
end

function f1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    return t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) + t0_t(p,q,kₐ,k₊,k₋,Eₐ,E₊))
end

function f1(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    return t1(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1(p,q,kₐ,k₊,k₋,Eₐ,E₊) + t0(p,q,kₐ,k₊,k₋,Eₐ,E₊))
end

function V̄_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)

    if p != q
        v_c = (((q/p)^ℓ-1)/((q/p)^(ℓ+1)-1)) * q * k₊⁺/(k₊⁺+k₋)
        else
            v_c = (ℓ)/(ℓ+1) * q * k₊⁺/(k₊⁺+k₋)
        end

    v_f = minimum([p,q])
    p0 = t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) + t0_t(p,q,kₐ,k₊,k₋,Eₐ,E₊))
    return v_f * (1-p0) + v_c * p0
end


function C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    return kₐ/(kₐ+k_d) * t1(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1(p,q,kₐ,k₊,k₋,Eₐ,E₊) + t0(p,q,kₐ,k₊,k₋,Eₐ,E₊))
end

function C_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    return kₐ/(kₐ+k_d) * t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) + t0_t(p,q,kₐ,k₊,k₋,Eₐ,E₊))
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
    if isinf(t0(p,q,kₐ,k₊,k₋,Eₐ,E₊))
        return 0
    else
        1 -  te(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1(p,q,kₐ,k₊,k₋,Eₐ,E₊)+t0(p,q,kₐ,k₊,k₋,Eₐ,E₊))
    end
end
function F_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ)
    if q > p
        p0 = (p/k₊) / ((q-p)/k₋ + p/k₊)
    else 
    p0 = 1
    end
    function te_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)
        if q > p
            p0 = (p/k₊) / ((q-p)/k₋ + p/k₊)
        else 
            p0 = 1
        end
        # p1 = p/(p+q) * ((kₐ+p)/(p+kₐ+q)) * (kₐ)/(kₐ + kₐ * exp(-Eₐ))
        if q*(k₊)/(k₊+k₋) > p
            return Inf
        else
            t0 = ((q/(p*k₋)+ (1/k₋ + 1/k₊) * (q-p)/p*k₊/k₋)* (exp(-k₋*(ℓₚ-ℓ)/(q-p))) +  (1-p0)*((q-p)/(p/k₊ - (q-p)/k₋)*1/k₋*(1/k₋ + 1/k₊)))
            return t0
        end
    end
    k₊⁺ = k₊ * exp(E₊)
    k_d = kₐ * exp(-Eₐ)
    k₊⁺ = k₊ * exp(E₊)
    if isinf(t0_t(p,q,kₐ,k₊,k₋,Eₐ,E₊))
        return 0
    else
        1 -  te_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)+t0_t(p,q,kₐ,k₊,k₋,Eₐ,E₊))
    end
end
function C₀(Eₐ)
    local kₐ = 1
    k_d = kₐ * exp(-Eₐ)
    return kₐ/(kₐ+k_d)
end
using CSV,DataFrames,Pipe
## testing for t0 and t1
begin
    q = 30
    kₐ = 100
    Eₐ = 3
    k₊ = 0.3
    k₋ = 0.4
    L = 335
    E₊ = 2
    ℓ = 4
    q̄ = q*k₊/(k₊+k₋)
    ps = [p for p in 2:0.2:100]
    t1s = [t1(p,q,kₐ,k₊,k₋,Eₐ,E₊) for p in ps]
    # ana_df = CSV.read("fig/analysis_df_analysis_p.csv",DataFrame)
    using Plots
    plotly()
    t0s = [t0_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) for p in ps]
    plot1=plot(ps,t0s, label = "t0_t")
    plot!(ps, [t0_em(p,q,kₐ,k₊,k₋,Eₐ,E₊) for p in ps], label = "t0_em")
    # plot!(ana_df.p[ana_df.p .≥ q̄],ana_df.T_0[ana_df.p .≥ q̄], label = "T_0 actual", line=:scatter)
    plot!(ps, [k₋^(-1) for p in ps],label = "k₋",line=:dash)
    plot!(ps, [q/(p*k₋) for p in ps], label= "q/(p*k₋)",line=:dash)
    # ylims!(0,15)
    yaxis!(:log)
    plot_2= plot(ps, [t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) for p in ps], label = "t1_t")
    plot!(ps, [t1_em(p,q,kₐ,k₊,k₋,Eₐ,E₊) for p in ps], label = "t1_em")
    # plot!(ana_df.p,ana_df.T_1, label = "T_1 actual",line=:scatter)
    plot!(ps, [ Tₙ(p,q,ℓ) * (1 + k₋/k₊*exp(E₊) + k₋*ℓ/p) for p in ps], label = "Tₙ")   
    yaxis!(:log)
    plot_3 = plot(ps, [t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t0_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)+t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊)) for p in ps], label="f1_t")
    plot!(ps, [t1_em(p,q,kₐ,k₊,k₋,Eₐ,E₊)/(t1_em(p,q,kₐ,k₊,k₋,Eₐ,E₊) +t0_em(p,q,kₐ,k₊,k₋,Eₐ,E₊)) for p in ps],label="f1_em")
    # # plot!(ana_df.p,ana_df.T_1./(ana_df.T_0.+ ana_df.T_1), label = "T_1/T_0 actual",line=:scatter)
    plot(plot1,plot_2,plot_3,layout=(3,1),size=(400,600))
    xlabel!("p")
    ylabel!("t0 and t1")
end
## testing for t0 and t1
begin
    Eₐs = [i for i in -2:0.25:7];
    p = 20;
    q = 30;
    E₊ = 2;
    kₐ = 100
    k₊ = 0.3
    k₋ = 0.4
    L = 335
    E₊ = 2
    ℓ = 4
    q̄ = q*k₊/(k₊+k₋)
    # ana_df = CSV.read("fig/analysis_df_analysis_p.csv",DataFrame)
    using Plots
    plotly()
    plot1=plot(Eₐs,[t0_em(p,q,kₐ,k₊,k₋,Eₐ,E₊) for Eₐ in Eₐs], label = "t0_em",line=:scatter)
    plot!(Eₐs, [t0_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) for Eₐ in Eₐs], label = "t0_t")
    # plot!(ana_df.p[ana_df.p .≥ q̄],ana_df.T_0[ana_df.p .≥ q̄], label = "T_0 actual", line=:scatter)
    plot!(Eₐs, [t1_em(p,q,kₐ,k₊,k₋,Eₐ,E₊) for Eₐ in Eₐs],label = "t1_em",line=:scatter)
    plot!(Eₐs, [t1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) for Eₐ in Eₐs], label = "t1_t")
    # ylims!(0,15)
    xlabel!("Eₐ")
    ylabel!("t0 and t1")
    Eₐs = [i for i in -2:0.25:7];
    p = 20;
    q = 30;
    E₊ = 2;
    Cₑₛₜ = [C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊) for Eₐ in Eₐs]
    plot2 = plot(Eₐs,Cₑₛₜ, label = "C₊")
    plot!(Eₐs, [C_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) for Eₐ in Eₐs], label = "C_t")
    df = CSV.read("fig/simu_df_Ebs_2.csv",DataFrame)
    plot!(df.Eᵦ, df.C, label="L=335", line=:scatter)
    xlabel!("Eₐ")
    ylabel!("C")

    Eₐs = [i for i in -2:0.25:7];
    p = 20;
    q = 30;
    E₊ = 2;
    ℓₚ = 27;
    Fₑₛₜ = [Fₜ₊(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ) for Eₐ in Eₐs]
    plot3=plot(Eₐs,Fₑₛₜ, label = "Fₜ₊")
    plot!(Eₐs, [F_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ) for Eₐ in Eₐs], label = "F_t")
    df = CSV.read("fig/simu_df_Ebs_2.csv",DataFrame)
    plot!(df.Eᵦ, df.fractions_T_protected, label="L=335", line=:scatter)
    xlabel!("Eₐ")
    ylabel!("F_T")
    df = CSV.read("fig/simu_df_Ebs_L_inf_2.csv",DataFrame)
    plot!(df.Eᵦ, df.fractions_T_protected, label="L→∞", line=:scatter)
    
    plot(plot1,plot2,plot3,layout=(3,1),size=(400,600))
end
## Testing about Coupling coefficient
begin
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
    p = [i for i in 2:0.2:100]
    Cₑₛₜ = [C₊(i,q,kₐ,k₊,k₋,Eₐ,E₊) for i in p]
    ana_df = CSV.read("fig/analysis_df_analysis_p.csv",DataFrame)

    plot1 = plot(p,Cₑₛₜ, label = "C₊")
    plot!(p, [C_t(i,q,kₐ,k₊,k₋,Eₐ,E₊) for i in p], label = "C_t")
    CSV.write("fig/plot_df_C_p.csv", DataFrame(p=p,eta=p./(q*k₊/(k₊+k₋)),C=[C_t(i,q,kₐ,k₊,k₋,Eₐ,E₊) for i in p]))

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

    k₊⁺ = k₊ * exp(E₊)
    # plot!(ana_df.p, ana_df.T_1 ./ (ana_df.T_0 .+ ana_df.T_1).*(1/(1 + k₋/k₊⁺ )*(kₐ)/(kₐ + kₐ * exp(-Eₐ))) , label="T_0/T_0+T_1", line=:scatter)
    plot!(p, [C₀(Eₐ) for i in p], label="C₀", line=:dash)
    xaxis!(:log)

#####################
    Eₐs = [i for i in -2:0.25:7];
    p = 20;
    q = 30;
    E₊ = 2;
    Cₑₛₜ = [C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊) for Eₐ in Eₐs]
    plot2 = plot(Eₐs,Cₑₛₜ, label = "C₊")
    plot!(Eₐs, [C_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) for Eₐ in Eₐs], label = "C_t")
    df = CSV.read("fig/simu_df_Ebs_2.csv",DataFrame)
    plot!(df.Eᵦ, df.C, label="L=335", line=:scatter)
    xlabel!("Eₐ")
    ylabel!("C")

    E₊s = [i for i in 0:0.1:5]
    p = 20;
    q = 30;
    Eᵦ = 3;
    Cₑₛₜ = [C₊(p,q,kₐ,k₊,k₋,Eₐ,E₊) for E₊ in E₊s]
    plot3 = plot(E₊s,Cₑₛₜ, label = "C₊")
    plot!(E₊s, [C_t(p,q,kₐ,k₊,k₋,Eₐ,E₊) for E₊ in E₊s], label = "C_t")
    df = CSV.read("fig/simu_df_E_cs.csv",DataFrame)
    plot!(df.E_c, df.C, label="L=335", line=:scatter)
    plot(plot1,plot2,plot3,layout=(3,1),size=(400,600))
    
end
## Testing about protection fraction
begin
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
    using Plots
    plotly()
    plot1 = plot(ps,Fₑₛₜ, label = "Fₜ₊")
    plot!(ps, [F_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ) for p in ps], label = "F_t")
    using DataFrames
    using CSV
    CSV.write("fig/plot_df_F_p.csv", DataFrame(p=ps,eta=ps./(q*k₊/(k₊+k₋)),F_T=[F_t(i,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ) for i in ps]))

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
    plot2=plot(Eₐs,Fₑₛₜ, label = "Fₜ₊")
    plot!(Eₐs, [F_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ) for Eₐ in Eₐs], label = "F_t")
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
    plot3 = plot(E₊s,Fₑₛₜ, label = "Fₜ₊")
    plot!(E₊s, [F_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ,ℓₚ) for E₊ in E₊s], label = "F_t")
    df = CSV.read("fig/simu_df_E_cs.csv",DataFrame)
    plot!(df.E_c, df.fractions_T_protected, label="L=335", line=:scatter)
    df = CSV.read("fig/simu_df_E_cs_L_inf.csv",DataFrame)

    xlabel!("E₊")
    ylabel!("F_T")
    plot(plot1,plot2,plot3,layout=(3,1),size=(400,600))
end
## Testing about the coupled velocity
begin
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
    plot1 = plot(ps,Vₑₛₜ, label = "V̄ₐ")
    plot!(ps, [V̄_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for p in ps], label = "V̄_t")

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
    plot2=plot(Eₐs,Vₑₛₜ, label = "V̄ₐ")
    plot!(Eₐs, [V̄_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for Eₐ in Eₐs], label = "V̄_t")
    df = CSV.read("fig/simu_df_Ebs_2.csv",DataFrame)
    plot!(df.Eᵦ, df.v_eff_translation, label="L=335", line=:scatter)
    xlabel!("Eₐ")
    ylabel!("C")

    E₊s = [i for i in 0:0.1:5]
    p = 20;
    q = 30;
    Eᵦ = 3;
    V̄ₐ(30,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ)
    Vₑₛₜ = [V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for E₊ in E₊s]
    plot3=plot(E₊s,Vₑₛₜ, label = "V̄ₐ")
    plot!(E₊s, [V̄_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for E₊ in E₊s], label = "V̄_t")
    df = CSV.read("fig/simu_df_E_cs.csv",DataFrame)
    plot!(df.E_c, df.v_eff_translation, label="L=335", line=:scatter)

    CSV.write("fig/plot_df_V_E_cs.csv", DataFrame(Ec=E₊s, V=[V̄_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for E₊ in E₊s]))

    E₊s = [i for i in 0:0.5:2]
    ps = [i for i in 13:0.2:60]
    Vₑₛₜ = [V̄ₐ(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for p in ps, E₊ in E₊s]
    plot4=plot(ps,Vₑₛₜ, label = "V̄ₐ")
    plot!(ps, [V̄_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for p in ps, E₊ in E₊s], label = "V̄_t")
    df = CSV.read("fig/merged_df_ps_var_E_cs.csv",DataFrame)
    for E₊ in E₊s
        plot!(df.p[df.Ec .== E₊], df.v_eff_transcriptions[df.Ec .== E₊], label="E₊=$E₊", line=:scatter)
    end
    xlabel!("p")
    plot5 = plot(ps,[f1(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for p in ps, E₊ in E₊s], label = "f1")
    plot!(ps,[f1_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ) for p in ps, E₊ in E₊s], label = "f1_t")
    plot_df = DataFrame(
        p = Float64[],
        V = Float64[],
        Ec = Float64[]
    )
    for E₊ in E₊s
        for p in ps
            push!(plot_df, [p, V̄_t(p,q,kₐ,k₊,k₋,Eₐ,E₊,ℓ), E₊])
        end
    end
    CSV.write("fig/plot_df_V_p_E_cs.csv", plot_df)
    plot(plot1,plot2,plot3,plot4,plot5,layout=(5,1),size=(400,1000))
end