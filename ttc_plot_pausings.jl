if isfile("simu_$(label).csv") & !overwrite    
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for K in qs
            cmd = `julia ttc_simu_base.jl $k $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $(k_stalling_0*K) $(k_unstalling_0*K) $k_ini_pausing $d $type $mode $N`
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
        s₀ = Int[], 
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

V_stalls = sort(unique(df.v_stalls))
# plot the effective velocity
sort!(df,[:v_stalls])
@pgf axis = SemiLogXAxis(
    {
        xlabel = "RNAP pausing rate \$k_-\$",
        ylabel = "effective velocity",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)
# group data by rate of translation
v_eff_transcription = Float64[]
std_eff_translation = Float64[]
v_eff_translation = Float64[]
std_eff_transcription = Float64[]
for V_stall in V_stalls
    T_transcriptions = df.T_transcription[df.v_stalls .== V_stall]
    T_translations = df.T_translation[df.v_stalls .== V_stall]
    push!(v_eff_transcription, L/mean(T_transcriptions))
    push!(std_eff_transcription, (L/(mean(T_transcriptions)-std(T_transcriptions)) - L/(mean(T_transcriptions)+std(T_transcriptions)) )/2)
    push!(v_eff_translation, L/mean(T_translations))
    push!(std_eff_translation, (L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2)
end
plot_line = Plot(Coordinates(V_stalls, v_eff_transcription; yerror= std_eff_transcription ))
legend_line = LegendEntry("transcription")
push!(axis, plot_line )
push!(axis, legend_line )
plot_line = Plot(Coordinates(V_stalls, v_eff_translation; yerror = std_eff_translation))
legend_line = LegendEntry("translation")
push!(axis, plot_line )
push!(axis, legend_line )

pgfsave("fig/effective_velocity_plot_$(label).tex",axis)
pgfsave("fig/effective_velocity_plot_$(label).svg",axis)

# plot the mean and std of the T_exposure
# sort!(df,[:v_translations])
@pgf axis = SemiLogXAxis(
    {
        xlabel = "RNAP pausing rate \$k_-\$",
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
for V_stall in V_stalls
    T_exposures = df.T_exposure[df.v_stalls .== V_stall]
    T_exposures_uncoupled = df.T_exposure_uncoupled[df.v_stalls .== V_stall]
    push!(mean_T_exposure, mean(T_exposures))
    push!(std_T_exposure, std(T_exposures))
    push!(mean_T_exposure_uncoupled, mean(T_exposures_uncoupled))
    push!(std_T_exposure_uncoupled, std(T_exposures_uncoupled))
end
plot_line = Plot(Coordinates(V_stalls, mean_T_exposure; yerror= std_T_exposure ))
legend_line = LegendEntry("exposure duration")
push!(axis, plot_line )
push!(axis, legend_line )
plot_line = Plot(Coordinates(V_stalls, mean_T_exposure_uncoupled; yerror = std_T_exposure_uncoupled))
legend_line = LegendEntry("exposure and uncoupled duration")
push!(axis, plot_line )
push!(axis, legend_line )

pgfsave("fig/exposure_duration_plot_$(label).tex",axis)
pgfsave("fig/exposure_duration_plot_$(label).svg",axis)

# plot the fraction of T_exposure in T_transcription
# sort!(df,[:v_translations])
@pgf axis = SemiLogXAxis(
    {
        xlabel = "RNAP pausing rate \$k_-\$",
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
for V_stall in V_stalls
    T_exposures = df.T_exposure[df.v_stalls .== V_stall]
    T_transcriptions = df.T_transcription[df.v_stalls .== V_stall]
    T_exposures_uncoupled = df.T_exposure_uncoupled[df.v_stalls .== V_stall]
    f_exposures =T_exposures ./ T_transcriptions
    f_exposures_uncoupled = T_exposures_uncoupled ./ T_transcriptions
    push!(fractions_T_exposure, mean(f_exposures))
    push!(std_fractions_T_exposure, std(f_exposures))
    push!(fractions_T_exposure_uncoupled, mean(f_exposures_uncoupled))
    push!(std_fractions_T_exposure_uncoupled, std(f_exposures_uncoupled))
end
plot_line = Plot(Coordinates(V_stalls, fractions_T_exposure; yerror= std_fractions_T_exposure ))
legend_line = LegendEntry("exposure fraction")
push!(axis, plot_line )
push!(axis, legend_line )
plot_line = Plot(Coordinates(V_stalls, fractions_T_exposure_uncoupled; yerror = std_fractions_T_exposure_uncoupled))
legend_line = LegendEntry("exposure and uncoupled fraction")
push!(axis, plot_line )
push!(axis, legend_line )

pgfsave("fig/exposure_fraction_plot_$(label).tex",axis)
pgfsave("fig/exposure_fraction_plot_$(label).svg",axis)


# plot the mean and std of the mean_distance
# sort!(df,[:v_translations])
@pgf axis = SemiLogXAxis(
    {
        xlabel = "RNAP pausing rate \$k_-\$",
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

for V_stall in V_stalls
    distances = df.mean_distance[df.v_stalls .== V_stall]
    push!(mean_distances, mean(distances))
    push!(std_distances, std(distances))
end
plot_line = Plot(Coordinates(V_stalls, mean_distances; yerror= std_distances ))
push!(axis, plot_line)

pgfsave("fig/mean_distance_plot_$(label).tex",axis)
pgfsave("fig/mean_distance_plot_$(label).svg",axis)
