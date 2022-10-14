# a schematic model for the LacZ complemention system
using CSV, DataFrames
T = 10 # delay
using NumericalIntegration
function ∫fdx(f,x,dx=1e-1)
    lin_space = collect(0:dx:x)
    if length(lin_space) > 1
    return integrate(lin_space,f.(lin_space))
    else
        return 0
    end
end
function ∫gdx(f,x,dx=1e-1)
    lin_space = collect(0:dx:x)
    if length(lin_space) > 1
    return sum(f.(lin_space))*dx
    else
        return 0
    end
end

X = collect(0:0.1:20)

function fit(X,GX)
    # GX = a*x + b
    a = (GX[end] - GX[end-5])/(X[end] - X[end-5])
    b = GX[end] - a*X[end]
    return a.*X .+ b
end

function intercept(X,GX)
    # GX = a*x + b
    a = (GX[end] - GX[end-5])/(X[end] - X[end-5])
    b = GX[end] - a*X[end]
    return -b/a
end

function taking_positive(X)
    return [maximum([0,x]) for x in X]
end

function model(f,X,label)
    ∫f(x) = ∫fdx(f,x)
    ∫∫f(x) = ∫fdx(∫f,x,1e-1)
    ∫∫∫f(x) = ∫fdx(∫∫f,x,1e-1)
    fX = f.(X)
    ∫fX = ∫f.(X)
    ∫∫fX = ∫∫f.(X)
    ∫∫∫fX = ∫∫∫f.(X)
    S_mRNA = ∫∫fX
    S_protein = sqrt.(∫∫∫fX)
    S_protein² = ∫∫∫fX
    fit_mRNA = fit(X,S_mRNA) |> taking_positive
    T_mRNA = intercept(X,S_mRNA)
    fit_protein = fit(X, S_protein) |> taking_positive
    T_protein = intercept(X,S_protein)
    temp_df = DataFrame(
        t = X,
        ρ = fX,
        S_mRNA = S_mRNA,
        S_protein = S_protein,
        S_protein² = S_protein²,
        fit_mRNA = fit_mRNA,
        T_mRNA = T_mRNA,
        fit_protein = fit_protein,
        T_protein = T_protein,
    )
    CSV.write("fig/sample_df_$(label).csv",temp_df)
    return (T_mRNA,T_protein)
end
f(x) = exp(-(x - T)^2/2)/√(2π)
h(x) = exp(-(x - T)^2)/√(π)
g(x) = (h(x) + h(x-5))/2
w₀(x) = (f(x) + f(x+3))/2
W₀ = ∫fdx(w₀,20)
w(x) = w₀(x)/W₀
Density = [
    f,
    g,
    w,
]
Label = [
    "unimodal",
    "distinct_bimodal",
    "adjacent_bimodal",
]
T_mRNA = []
T_protein = []
for (f,label) in zip(Density,Label)
    t_mRNA, t_protein=model(f,X,label)
    push!(T_mRNA,t_mRNA)
    push!(T_protein,t_protein)
end
df = DataFrame(
    Label = Label,
    T_mRNA = T_mRNA,
    T_protein = T_protein,
)
CSV.write("fig/sample_df_Ts.csv",df)

