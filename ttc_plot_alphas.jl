if isfile("simu_$(label).csv") & !overwrite
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for α in αs
        cmd = `julia ttc_simu_base.jl $k $p $L $ℓ $(α) $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $mode $N`
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
        );

    for f in readdir("./data/",join=true)
        temp_df = CSV.read(f,DataFrame)
        for i in 1:length(temp_df[:,1])
            push!(df, temp_df[i,:])
        end
        rm(f)
    end
    CSV.write("simu_$(label).csv",df)
end
using LaTeXStrings

if @isdefined approx_flag
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
    for α in αs
        T_transcriptions = df.T_transcription[df.α .== α]
        T_translations = df.T_translation[df.α .== α]
        push!(v_eff_transcription, L/mean(T_transcriptions))
        push!(std_eff_transcription, (L/(mean(T_transcriptions)-std(T_transcriptions)) - L/(mean(T_transcriptions)+std(T_transcriptions)) )/2)
        push!(v_eff_translation, L/mean(T_translations))
        push!(std_eff_translation, (L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2)
        T_protecteds = df.T_exposure[df.α .== α]
        T_transcriptions = df.T_transcription[df.α .== α]
        T_protecteds_uncoupled = df.T_exposure_uncoupled[df.α .== α]
        f_protecteds =T_protecteds ./ T_transcriptions
        f_protecteds_uncoupled = T_protecteds_uncoupled ./ T_transcriptions
        push!(fractions_T_protected, 1-mean(f_protecteds))
        push!(std_fractions_T_protected, std(f_protecteds))
        push!(fractions_T_protected_uncoupled, mean(f_protecteds_uncoupled))
        push!(std_fractions_T_protected_uncoupled, std(f_protecteds_uncoupled))
        push!(v_eff_est, V̄(p,k,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c,ℓ))
        push!(F_T_est, 𝔼Fₜ(p,k,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c,ℓ,27))
        push!(C₊_est, C₊(p,k,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c))
        push!(Cₐ_est, Cₐ(p,k,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c, α, L))
        push!(C, mean(df.s[df.α .== α]))
    end

    plot_df = DataFrame(
            α=αs, 
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


# plot the effective velocity
sort!(df,[:α])
@pgf axis = SemiLogXAxis(
    {
        xlabel = "translation initiation rate \$\\alpha\$",
        ylabel = "effective velocity",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        # legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)
# group data by rate of translation
v_eff_transcription = Float64[]
std_eff_translation = Float64[]
v_eff_translation = Float64[]
std_eff_transcription = Float64[]
for α in αs
    T_transcriptions = df.T_transcription[df.α .== α]
    T_translations = df.T_translation[df.α .== α]
    push!(v_eff_transcription, L/mean(T_transcriptions))
    push!(std_eff_transcription, (L/(mean(T_transcriptions)-std(T_transcriptions)) - L/(mean(T_transcriptions)+std(T_transcriptions)) )/2)
    push!(v_eff_translation, L/mean(T_translations))
    push!(std_eff_translation, (L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2)
end
plot_line = Plot(Coordinates(αs, v_eff_transcription; yerror= std_eff_transcription ))
legend_line = LegendEntry("transcription")
push!(axis, plot_line )
push!(axis, legend_line )
plot_line = Plot(Coordinates(αs, v_eff_translation; yerror = std_eff_translation))
legend_line = LegendEntry("translation")
push!(axis, plot_line )
push!(axis, legend_line )
effective_velocity_df = DataFrame(
    alpha = αs,
    v_eff_transcription = v_eff_transcription,
    std_eff_transcription = std_eff_transcription,
    v_eff_translation = v_eff_translation,
    std_eff_translation = std_eff_translation
)
CSV.write("fig/effective_velocity_df_$(label).csv",effective_velocity_df)
pgfsave("fig/effective_velocity_plot_$(label).tex",axis)
pgfsave("fig/effective_velocity_plot_$(label).svg",axis)

# plot the mean and std of the T_exposure
# sort!(df,[:v_translations])
@pgf axis = SemiLogXAxis(
    {
        xlabel = "translation initiation rate \$\\alpha\$",
        ylabel = "duration / s",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)
# group data by rate of translation
mean_T_exposure = Float64[]
std_T_exposure = Float64[]
mean_T_exposure_uncoupled = Float64[]
std_T_exposure_uncoupled = Float64[]
for α in αs
    T_exposures = df.T_exposure[df.α .== α]
    T_exposures_uncoupled = df.T_exposure_uncoupled[df.α .== α]
    push!(mean_T_exposure, mean(T_exposures))
    push!(std_T_exposure, std(T_exposures))
    push!(mean_T_exposure_uncoupled, mean(T_exposures_uncoupled))
    push!(std_T_exposure_uncoupled, std(T_exposures_uncoupled))
end
plot_line = Plot(Coordinates(αs, mean_T_exposure; yerror= std_T_exposure ))
legend_line = LegendEntry("exposure duration")
push!(axis, plot_line )
push!(axis, legend_line )
plot_line = Plot(Coordinates(αs, mean_T_exposure_uncoupled; yerror = std_T_exposure_uncoupled))
legend_line = LegendEntry("exposure and uncoupled duration")
push!(axis, plot_line )
push!(axis, legend_line )



pgfsave("fig/exposure_duration_plot_$(label).tex",axis)
pgfsave("fig/exposure_duration_plot_$(label).svg",axis)

# plot the fraction of T_exposure in T_transcription
# sort!(df,[:v_translations])
@pgf axis = SemiLogXAxis(
    {
        xlabel = "translation initiation rate \$\\alpha\$",
        ylabel = "fraction",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)
# group data by rate of translation
fractions_T_exposure = Float64[]
fractions_T_exposure_uncoupled = Float64[]
std_fractions_T_exposure = Float64[]
std_fractions_T_exposure_uncoupled = Float64[]
for α in αs
    T_exposures = df.T_exposure[df.α .== α]
    T_transcriptions = df.T_transcription[df.α .== α]
    T_exposures_uncoupled = df.T_exposure_uncoupled[df.α .== α]
    f_exposures =T_exposures ./ T_transcriptions
    f_exposures_uncoupled = T_exposures_uncoupled ./ T_transcriptions
    push!(fractions_T_exposure, mean(f_exposures))
    push!(std_fractions_T_exposure, std(f_exposures))
    push!(fractions_T_exposure_uncoupled, mean(f_exposures_uncoupled))
    push!(std_fractions_T_exposure_uncoupled, std(f_exposures_uncoupled))
end
plot_line = Plot(Coordinates(αs, fractions_T_exposure; yerror= std_fractions_T_exposure ))
legend_line = LegendEntry("exposure fraction")
push!(axis, plot_line )
push!(axis, legend_line )
plot_line = Plot(Coordinates(αs, fractions_T_exposure_uncoupled; yerror = std_fractions_T_exposure_uncoupled))
legend_line = LegendEntry("exposure and uncoupled fraction")
push!(axis, plot_line )
push!(axis, legend_line )
protected_fraction_df = DataFrame(
    alpha = αs,
    fractions_T_protected = - fractions_T_exposure .+ 1,
    std_fractions_T_protected = std_fractions_T_exposure
)
CSV.write("fig/protected_fraction_df_$(label).csv",protected_fraction_df)
pgfsave("fig/exposure_fraction_plot_$(label).tex",axis)
pgfsave("fig/exposure_fraction_plot_$(label).svg",axis)

# plot the mean and std of the mean_distance
# sort!(df,[:v_translations])
@pgf axis = SemiLogXAxis(
    {
        xlabel = "translation initiation rate \$\\alpha\$",
        ylabel = "distance / codon",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)
# group data by rate of translation
mean_distances = Float64[]
std_distances = Float64[]
for α in αs
    distances = df.mean_distance[df.α .== α]
    push!(mean_distances, mean(distances))
    push!(std_distances, std(distances))
end
plot_line = Plot(Coordinates(αs, mean_distances; yerror= std_distances ))
# legend_line = LegendEntry("distance")
push!(axis, plot_line )
# push!(axis, legend_line )

pgfsave("fig/mean_distance_plot_$(label).tex",axis)
pgfsave("fig/mean_distance_plot_$(label).svg",axis)

# merged_df


merged_df = DataFrame(
    alpha = Float64[],
    v_translations = Float64[],
    v_transcriptions = Float64[],
    fraction_protected = Float64[],
    std_fraction_protected = Float64[],
    v_eff_translations = Float64[],
    v_eff_transcriptions = Float64[],
    std_eff_translations = Float64[],
    std_eff_transcriptions = Float64[],
    efficiency = Float64[],
    std_efficiency = Float64[],
)

for α in αs
    v_eff_transcriptions = (L)./df.T_transcription[df.α .== α]
    v_eff_translations = (L)./df.T_translation[df.α .== α]
    T_exposures = df.T_exposure[df.α .== α]
    T_transcriptions = df.T_transcription[df.α .== α]
    f_exposures =T_exposures ./ T_transcriptions
    fraction_T_protected = - f_exposures .+ 1
    push!(merged_df,
        [
            α,
            p,
            k,
            mean(fraction_T_protected),
            std(fraction_T_protected),
            mean(v_eff_translations),
            mean(v_eff_transcriptions),
            std(v_eff_translations),
            std(v_eff_transcriptions),
            mean(v_eff_translations) ./ p,
            std(v_eff_translations) ./ p,
        ]
    )
end

CSV.write(
    "fig/merged_df_$(label).csv",
    merged_df
    )

@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        ylabel = "ribosome pushing efficiency \$ \\eta = \\bar V_{\\rm rib}/p\$",
        xlabel = "mean fraction of protected duration",
        # "error bars/y dir=both",
        # "error bars/y explicit",
        # "error bars/x dir=both",
        # "error bars/x explicit",
    },
)

plot_line = Plot(Coordinates(merged_df.fraction_protected, merged_df.efficiency))
# legend_line = LegendEntry("efficiency")
push!(axis,plot_line)

pgfsave("fig/merged_plot_$(label).tex",axis)
pgfsave("fig/merged_plot_$(label).svg",axis)