

if isfile("simu_$(label).csv") & append_flag
    df = CSV.read("simu_$(label).csv",DataFrame)
    cmds = []
    for p in ps
        cmd = `julia ttc_simu_base.jl $k $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $mode $N`
        push!(cmds, cmd)
    end
    @showprogress 1 pmap(run, cmds)
    for f in readdir("./data/",join=true)
        temp_df = CSV.read(f,DataFrame)
        for i in 1:length(temp_df[:,1])
            push!(df, temp_df[i,:])
        end
        rm(f)
    end
    CSV.write("simu_$(label).csv",df)
elseif isfile("simu_$(label).csv") & !overwrite    
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        cmd = `julia ttc_simu_base.jl $k $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $mode $N`
        push!(cmds, cmd)
    end
    @showprogress 1 pmap(run, cmds)

    df = DataFrame(
        k_couple = Float64[], 
        k_uncouple = Float64[], 
        v_translations = Float64[], 
        v_transcriptions = Float64[], 
        v_stalls = Float64[], 
        v_unstalls = Float64[], 
        α = Float64[],
        k_ini_pausings = Float64[] , 
        L = [], 
        ℓ = [], 
        Eᵦ = Float64[], 
        E_c = Float64[], 
        x₀ = Int[], 
        y₀ = Int[], 
        s = Int[], 
        p₀ = Int[], 
        type = String[], 
        T_transcription = Float64[],
        T_translation = Float64[], 
        T_exposure = Float64[], 
        T_exposure_uncoupled = Float64[], 
        d = Float64[], 
        mean_distance = Float64[],
        mode = String[],
        );    # end

    for f in readdir("./data/",join=true)
        temp_df = CSV.read(f,DataFrame)
        for i in 1:length(temp_df[:,1])
            push!(df, temp_df[i,:])
        end
        rm(f)
    end
    CSV.write("simu_$(label).csv",df)
end

if @isdefined approx_flag
    if approx_flag == false
        v_eff_transcription = Float64[]
        std_eff_translation = Float64[]
        v_eff_translation = Float64[]
        std_eff_transcription = Float64[]
        fractions_T_protected = Float64[]
        fractions_T_protected_uncoupled = Float64[]
        std_fractions_T_protected = Float64[]
        std_fractions_T_protected_uncoupled = Float64[]
        v_eff_est = Float64[]
        F_T_est = Float64[]
        C = Float64[]
        C₊_est = Float64[]
        Cₐ_est = Float64[]
        for p in ps
            T_transcriptions = df.T_transcription[df.v_translations .== p]
            T_translations = df.T_translation[df.v_translations .== p]
            push!(v_eff_transcription, L/mean(T_transcriptions))
            push!(std_eff_transcription, (L/(mean(T_transcriptions)-std(T_transcriptions)) - L/(mean(T_transcriptions)+std(T_transcriptions)) )/2)
            push!(v_eff_translation, L/mean(T_translations))
            push!(std_eff_translation, (L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2)
            T_protecteds = df.T_exposure[df.v_translations .== p]
            T_transcriptions = df.T_transcription[df.v_translations .== p]
            T_protecteds_uncoupled = df.T_exposure_uncoupled[df.v_translations .== p]
            f_protecteds =T_protecteds ./ T_transcriptions
            f_protecteds_uncoupled = T_protecteds_uncoupled ./ T_transcriptions
            push!(fractions_T_protected, 1-mean(f_protecteds))
            push!(std_fractions_T_protected, std(f_protecteds))
            push!(fractions_T_protected_uncoupled, mean(f_protecteds_uncoupled))
            push!(std_fractions_T_protected_uncoupled, std(f_protecteds_uncoupled))
        end
        
        plot_df = DataFrame(
                p=ps, 
                v_eff_transcription=v_eff_transcription, 
                std_eff_transcription=std_eff_transcription, 
                v_eff_translation=v_eff_translation, 
                std_eff_translation=std_eff_translation, 
                fractions_T_protected=fractions_T_protected,
                std_fractions_T_protected=std_fractions_T_protected,
                fractions_T_protected_uncoupled=fractions_T_protected_uncoupled,
                std_fractions_T_protected_uncoupled=std_fractions_T_protected_uncoupled,
            )
        CSV.write("fig/simu_df_$(label).csv",plot_df)
    else
        include("ttc_approx.jl")
        # group data by rate of translation
        v_eff_transcription = Float64[]
        std_eff_translation = Float64[]
        v_eff_translation = Float64[]
        std_eff_transcription = Float64[]
        fractions_T_protected = Float64[]
        fractions_T_protected_uncoupled = Float64[]
        std_fractions_T_protected = Float64[]
        std_fractions_T_protected_uncoupled = Float64[]
        v_eff_est = Float64[]
        F_T_est = Float64[]
        C = Float64[]
        C₊_est = Float64[]
        Cₐ_est = Float64[]
        for p in ps
            S = df.s[df.v_translations .== p]
            C_simu = mean(S)
            T_transcriptions = df.T_transcription[df.v_translations .== p]
            T_translations = df.T_translation[df.v_translations .== p]
            push!(v_eff_transcription, L/mean(T_transcriptions))
            push!(std_eff_transcription, (L/(mean(T_transcriptions)-std(T_transcriptions)) - L/(mean(T_transcriptions)+std(T_transcriptions)) )/2)
            push!(v_eff_translation, L/mean(T_translations))
            push!(std_eff_translation, (L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2)
            T_protecteds = df.T_exposure[df.v_translations .== p]
            T_transcriptions = df.T_transcription[df.v_translations .== p]
            T_protecteds_uncoupled = df.T_exposure_uncoupled[df.v_translations .== p]
            f_protecteds =T_protecteds ./ T_transcriptions
            f_protecteds_uncoupled = T_protecteds_uncoupled ./ T_transcriptions
            push!(fractions_T_protected, 1-mean(f_protecteds))
            push!(std_fractions_T_protected, std(f_protecteds))
            push!(fractions_T_protected_uncoupled, mean(f_protecteds_uncoupled))
            push!(std_fractions_T_protected_uncoupled, std(f_protecteds_uncoupled))
            push!(v_eff_est, V̄ₐ(p,k,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c,ℓ))
            push!(F_T_est, Fₜ₊(p,k,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c,ℓ,27))
            push!(C₊_est, C₊(p,k,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c))
            push!(Cₐ_est, C₊(p,k,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c))
            push!(C, C_simu)
        end

        plot_df = DataFrame(
                p=ps, 
                v_eff_transcription=v_eff_transcription, 
                std_eff_transcription=std_eff_transcription, 
                v_eff_translation=v_eff_translation, 
                std_eff_translation=std_eff_translation, 
                fractions_T_protected=fractions_T_protected,
                std_fractions_T_protected=std_fractions_T_protected,
                fractions_T_protected_uncoupled=fractions_T_protected_uncoupled,
                std_fractions_T_protected_uncoupled=std_fractions_T_protected_uncoupled,
                v_eff_est=v_eff_est,
                F_T_est=F_T_est,
                C₊_est = C₊_est,
                Cₐ_est = Cₐ_est,
                C = C,
            )
        CSV.write("fig/simu_df_$(label).csv",plot_df)
    end   
else
    v_eff_transcription = Float64[]
    std_eff_translation = Float64[]
    v_eff_translation = Float64[]
    std_eff_transcription = Float64[]
    fractions_T_protected = Float64[]
    fractions_T_protected_uncoupled = Float64[]
    std_fractions_T_protected = Float64[]
    std_fractions_T_protected_uncoupled = Float64[]
    v_eff_est = Float64[]
    F_T_est = Float64[]
    C = Float64[]
    C₊_est = Float64[]
    Cₐ_est = Float64[]
    for p in ps
        T_transcriptions = df.T_transcription[df.v_translations .== p]
        T_translations = df.T_translation[df.v_translations .== p]
        push!(v_eff_transcription, L/mean(T_transcriptions))
        push!(std_eff_transcription, (L/(mean(T_transcriptions)-std(T_transcriptions)) - L/(mean(T_transcriptions)+std(T_transcriptions)) )/2)
        push!(v_eff_translation, L/mean(T_translations))
        push!(std_eff_translation, (L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2)
        T_protecteds = df.T_exposure[df.v_translations .== p]
        T_transcriptions = df.T_transcription[df.v_translations .== p]
        T_protecteds_uncoupled = df.T_exposure_uncoupled[df.v_translations .== p]
        f_protecteds =T_protecteds ./ T_transcriptions
        f_protecteds_uncoupled = T_protecteds_uncoupled ./ T_transcriptions
        push!(fractions_T_protected, 1-mean(f_protecteds))
        push!(std_fractions_T_protected, std(f_protecteds))
        push!(fractions_T_protected_uncoupled, mean(f_protecteds_uncoupled))
        push!(std_fractions_T_protected_uncoupled, std(f_protecteds_uncoupled))
    end
    
    plot_df = DataFrame(
            p=ps, 
            v_eff_transcription=v_eff_transcription, 
            std_eff_transcription=std_eff_transcription, 
            v_eff_translation=v_eff_translation, 
            std_eff_translation=std_eff_translation, 
            fractions_T_protected=fractions_T_protected,
            std_fractions_T_protected=std_fractions_T_protected,
            fractions_T_protected_uncoupled=fractions_T_protected_uncoupled,
            std_fractions_T_protected_uncoupled=std_fractions_T_protected_uncoupled,
        )
    CSV.write("fig/simu_df_$(label).csv",plot_df)
end


using StatsBase, LinearAlgebra
using Pipe

function pdf_corr(X,Y)
    ∑(x) = sum(x)
    min = @pipe minimum([X;Y]) |> floor(Int,_*10)/10
    max = @pipe maximum([X;Y]) |> ceil(Int,_*10)/10
    dx = (max - min)/minimum([length(X),length(Y)])*√(minimum([length(X),length(Y)]))
    h_x = @pipe fit(Histogram, X, min:dx:max) |> 
        normalize(_, mode=:pdf)
    h_y = @pipe fit(Histogram, Y, min:dx:max) |> 
        normalize(_, mode=:pdf)
    f_x = h_x.weights
    f_y = h_y.weights
    ∫fgdx(f,g,dx) = ∑(f.*g)*dx
    corr = ∫fgdx(f_x,f_y,dx)/√(∫fgdx(f_x,f_x,dx)*∫fgdx(f_y,f_y,dx))
    return corr
end

# group df by the value of ps and qs
gdf = @pipe df |> 
    rename(_, :v_translations => :p) |> rename(_, :v_transcriptions => :q) |>
    sort(_, [:p, :q]) |> 
    groupby(_, [:p, :q])

function empirical_distribution(df)
    T_transcriptions = df.T_transcription 
    T_translations = df.T_translation
    tmin = 0.0
    tmax = 70.0
    dt = 0.5
    function renormalize(x)
        if maximum(x) == 0
            return x
        else
            return x./maximum(x)
        end
    end
    ρ_transcription = @pipe fit(Histogram, T_transcriptions, tmin:dt:tmax) |> 
        normalize(_, mode=:pdf) |> replace(_.weights, NaN => 0) |> renormalize
    ρ_translation = @pipe fit(Histogram, T_translations, tmin:dt:tmax) |>
        normalize(_, mode=:pdf) |> replace(_.weights, NaN => 0) |> renormalize
    tmax = 60.0
    local p = df.p[1]
    local q = df.q[1]
    edf = DataFrame(
        t = collect(tmin:dt:tmax-dt),
        pdf_transcription = ρ_transcription[1:length(collect(tmin:dt:tmax-dt))],
        pdf_translation = ρ_translation[1:length(collect(tmin:dt:tmax-dt))],
        p = [p for i in tmin:dt:tmax-dt],
        q = [q for i in tmin:dt:tmax-dt],
    )
    CSV.write("data/pdf_$(p)_$(q).csv",edf)
    # return ρ_transcription, ρ_translation
end

for i in 1:length(gdf)
    empirical_distribution(gdf[i])
end

pdf = DataFrame(
    t = Float64[],
    pdf_transcription = Float64[],
    pdf_translation = Float64[],
    p = Float64[],
    q = Float64[],
)
for f in readdir("./data/",join=true)
    temp_df = CSV.read(f,DataFrame)
    for i in 1:length(temp_df[:,1])
        push!(pdf, temp_df[i,:])
    end
    rm(f)
end
sort!(pdf,[:p,:q])
CSV.write("fig/pdf_$(label).csv",pdf)

using LaTeXStrings
# plot the effective velocity
sort!(df,[:v_translations])
@pgf axis = Axis(
    {   
        width = "3in",
        height = "3in",
        xlabel = "ribosome translocation rate \$p\$",
        ylabel = "effective velocity {\\color{red}\$\\overline V_{\\rm RNAP}\$}, {\\color{blue}\$\\overline V_{\\rm rib}\$}, \$\\overline V_{\\rm est}\$",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "south east",
    },
    VLine({ dotted, red }, k*(k_unstalling_0)/(k_stalling_0+k_unstalling_0))
)
t = @pgf Table({
    x = "p",
    y = "v_eff_transcription",
    col_sep = "comma"},"simu_df_$(label).csv")
print_tex(t)
plot_line = @pgf Plot({color="red",mark="o",only_marks},t)
legend_line = LegendEntry("transcription")
push!(axis, plot_line )
push!(axis, legend_line )
t = @pgf Table({
    x = "p",
    y = "v_eff_translation",
    col_sep = "comma"},"simu_df_$(label).csv")
plot_line = @pgf Plot({color="blue",mark="square",only_marks},t)
legend_line = LegendEntry("translation")
push!(axis, plot_line )
push!(axis, legend_line )
t = @pgf Table({
    x = "p",
    y = "v_eff_est",
    col_sep = "comma"},"simu_df_$(label).csv")
plot_line = @pgf Plot({color="black",no_marks,restrict_x_to_domain="$(k*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)):$(maximum(ps))"},t)
legend_line = LegendEntry("est")
push!(axis, plot_line )
push!(axis, legend_line )

pgfsave("fig/effective_velocity_plot_$(label).tex",axis)
# pgfsave("fig/effective_velocity_plot_$(label).svg",axis)

# plot the fraction of T_exposure in T_transcription
# sort!(df,[:v_translations])
@pgf axis = Axis(
    {
        xlabel = "ribosome translocation rate \$p\$",
        ylabel = "protected fraction \$ F_T\$",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)
t = @pgf Table({
    x = "p",
    y = "fractions_T_protected",
    col_sep = "comma"},"simu_df_$(label).csv")

plot_line = @pgf Plot({color="blue",only_marks},t)
push!(axis, plot_line )
t = @pgf Table({
    x = "p",
    y = "F_T_est",
    col_sep = "comma"},"simu_df_$(label).csv")

plot_line = @pgf PlotInc({restrict_x_to_domain="$(k*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)):$(maximum(ps))"},t)
push!(axis, plot_line )
push!(axis, legend_line )

pgfsave("fig/exposure_fraction_plot_$(label).tex",axis)

