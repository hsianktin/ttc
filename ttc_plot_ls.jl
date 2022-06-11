if isfile("simu_$(label).csv") & !overwrite    
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for ℓ in ℓs
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

# plot the effective velocity
sort!(df,[:ℓ])
@pgf axis = Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "interaction length \$\\ell\$",
        ylabel = "effective velocity",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north east"
    },
    # VLine({ dotted, red }, k)

    )
# group data by rate of translation
v_eff_transcription = Float64[]
std_eff_translation = Float64[]
v_eff_translation = Float64[]
std_eff_transcription = Float64[]
for ℓ in ℓs
    v_eff_transcriptions = (L)./df.T_transcription[df.ℓ .== ℓ]
    v_eff_translations = (L)./df.T_translation[df.ℓ .== ℓ]
    push!(v_eff_transcription, mean(v_eff_transcriptions))
    push!(std_eff_translation, std(v_eff_translations))
    push!(v_eff_translation, mean(v_eff_translations))
    push!(std_eff_transcription, std(v_eff_transcriptions))
end
plot_df = DataFrame(
    ℓs = ℓs, 
    v_eff_translation = v_eff_translation, 
    std_eff_translation = std_eff_translation, 
    mse_eff_translation = std_eff_translation/sqrt(1000),
    v_eff_transcription = v_eff_transcription, 
    std_eff_transcription = std_eff_transcription,
    mse_eff_transcription = std_eff_transcription/sqrt(1000),
    )
CSV.write("fig/effective_velocity_df_$(label).csv",plot_df)
plot_df = DataFrame(
    ℓs = ℓs[[i for i in 1:3:100]], 
    v_eff_translation = v_eff_translation[[i for i in 1:3:100]], 
    std_eff_translation = std_eff_translation[[i for i in 1:3:100]], 
    mse_eff_translation = std_eff_translation[[i for i in 1:3:100]]/sqrt(1000),
    v_eff_transcription = v_eff_transcription[[i for i in 1:3:100]], 
    std_eff_transcription = std_eff_transcription[[i for i in 1:3:100]],
    mse_eff_transcription = std_eff_transcription[[i for i in 1:3:100]]/sqrt(1000),
    )
CSV.write("fig/effective_velocity_df_$(label)_one_third.csv",plot_df)


t = @pgf Table({x = "ℓs", y = "v_eff_transcription", "y error"="mse_eff_transcription", "col sep"="comma"}, "effective_velocity_df_$(label).csv")
plot_line = Plot(t)
legend_line = LegendEntry("transcription")
push!(axis, plot_line )
push!(axis, legend_line )
t = @pgf Table({x = "ℓs", y = "v_eff_translation", "y error"="mse_eff_translation", "col sep"="comma"}, "effective_velocity_df_$(label).csv")
plot_line = Plot(t)
legend_line = LegendEntry("translation")
push!(axis, plot_line )
push!(axis, legend_line )

pgfsave("fig/effective_velocity_plot_$(label).tex",axis)

## alternative 

# group data by rate of translation
v_eff_transcription = Float64[]
std_eff_translation = Float64[]
v_eff_translation = Float64[]
std_eff_transcription = Float64[]
for ℓ in ℓs
    v_eff_transcriptions = (L)./df.T_transcription[df.ℓ .== ℓ]
    v_eff_translations = (L)./df.T_translation[df.ℓ .== ℓ]
    push!(v_eff_transcription, mean(v_eff_transcriptions))
    push!(std_eff_translation, std(v_eff_translations))
    push!(v_eff_translation, mean(v_eff_translations))
    push!(std_eff_transcription, std(v_eff_transcriptions))
end
plot_df = DataFrame(
    ℓs = ℓs,
    ℓs_shift = ℓs.+ 0.5,
    v_eff_translation = v_eff_translation, 
    std_eff_translation = std_eff_translation, 
    mse_eff_translation = std_eff_translation/sqrt(10000),
    v_eff_transcription = v_eff_transcription, 
    std_eff_transcription = std_eff_transcription,
    mse_eff_transcription = std_eff_transcription/sqrt(10000),
    )
CSV.write("fig/effective_velocity_df_$(label).csv",plot_df)

plot_df = DataFrame(
    ℓs = ℓs[[i for i in 1:3:100]], 
    ℓs_shift = ℓs[[i for i in 1:3:100]].+ 0.5,
    v_eff_translation = v_eff_translation[[i for i in 1:3:100]], 
    std_eff_translation = std_eff_translation[[i for i in 1:3:100]], 
    mse_eff_translation = std_eff_translation[[i for i in 1:3:100]]/sqrt(10000),
    v_eff_transcription = v_eff_transcription[[i for i in 1:3:100]], 
    std_eff_transcription = std_eff_transcription[[i for i in 1:3:100]],
    mse_eff_transcription = std_eff_transcription[[i for i in 1:3:100]]/sqrt(10000),
    )
CSV.write("fig/effective_velocity_df_$(label)_one_third.csv",plot_df)
t1 = @pgf Table({x = "ℓs", y = "v_eff_transcription", "y error"="std_eff_transcription", "col sep"="comma"}, "effective_velocity_df_$(label)_one_third.csv")
t2 = @pgf Table({x = "ℓs_shift", y = "v_eff_translation", "y error"="std_eff_translation", "col sep"="comma"}, "effective_velocity_df_$(label)_one_third.csv")
t3 = @pgf Table({x = "ℓs", y = "v_eff_transcription", "col sep"="comma"}, "effective_velocity_df_$(label).csv")
t4 = @pgf Table({x = "ℓs_shift", y = "v_eff_translation", "col sep"="comma"}, "effective_velocity_df_$(label).csv")

@pgf axis = Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "interaction length \$\\ell\$",
        ylabel = "effective velocity",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north east",
        ymin = 12,
        ymax = 26,
    },
    # VLine({ dotted, red }, k)
    Plot({color = "blue",mark = "*"},t3),
    LegendEntry("transcription"),
    Plot({color = "red",mark = "*"},t4),
    LegendEntry("transcription"),
    Plot(
        {"only_marks",color = "blue",mark = "*"},
        t1),
    Plot(
        {"only_marks",color = "red",mark = "*"},
        t2),
    )

# plot_line = Plot({"only_marks",color = "blue"},t1)
# legend_line = LegendEntry("transcription")
# push!(axis, plot_line )
# push!(axis, legend_line )
# plot_line = Plot({"only_marks",color = "red"},t2)
# legend_line = LegendEntry("translation")
# push!(axis, plot_line )
# push!(axis, legend_line )
# plot_line = Plot({color = "blue"},t3)
# push!(axis, plot_line )
# plot_line = Plot({color = "red"},t4)
# push!(axis, plot_line )


pgfsave("fig/effective_velocity_plot_$(label)_alt.tex",axis)


# pgfsave("fig/effective_velocity_plot_$(label).svg",axis)

# plot the mean and std of the T_exposure
# sort!(df,[:v_translations])
# @pgf axis = Axis(
#     {
#         xlabel = "interaction length \$\\ell\$",
#         ylabel = "duration / s",
#         grid = "major",
#         "error bars/y dir=both",
#         "error bars/y explicit",
#         legend_pos  = "north west"
#     },
#     # VLine({ dotted, red }, k)
# )
# # group data by rate of translation
# mean_T_exposure = Float64[]
# std_T_exposure = Float64[]
# mean_T_exposure_uncoupled = Float64[]
# std_T_exposure_uncoupled = Float64[]
# for ℓ in ℓs
#     T_exposures = df.T_exposure[df.ℓ .== ℓ]
#     T_exposures_uncoupled = df.T_exposure_uncoupled[df.ℓ .== ℓ]
#     push!(mean_T_exposure, mean(T_exposures))
#     push!(std_T_exposure, std(T_exposures))
#     push!(mean_T_exposure_uncoupled, mean(T_exposures_uncoupled))
#     push!(std_T_exposure_uncoupled, std(T_exposures_uncoupled))
# end
# plot_line = Plot(Coordinates(ℓs, mean_T_exposure; yerror= std_T_exposure ))
# legend_line = LegendEntry("exposure duration")
# push!(axis, plot_line )
# push!(axis, legend_line )
# plot_line = Plot(Coordinates(ℓs, mean_T_exposure_uncoupled; yerror = std_T_exposure_uncoupled))
# legend_line = LegendEntry("exposure and uncoupled duration")
# push!(axis, plot_line )
# push!(axis, legend_line )

# pgfsave("fig/exposure_duration_plot_$(label).tex",axis)
# pgfsave("fig/exposure_duration_plot_$(label).svg",axis)

# plot the fraction of T_exposure in T_transcription
# sort!(df,[:v_translations])
@pgf axis = Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "interaction length \$\\ell\$",
        ylabel = "fraction",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north west"
    },
    )
# group data by rate of translation
fractions_T_exposure = Float64[]
fractions_T_exposure_uncoupled = Float64[]
std_fractions_T_exposure = Float64[]
std_fractions_T_exposure_uncoupled = Float64[]
for ℓ in ℓs
    T_exposures = df.T_exposure[df.ℓ .== ℓ]
    T_transcriptions = df.T_transcription[df.ℓ .== ℓ]
    T_exposures_uncoupled = df.T_exposure_uncoupled[df.ℓ .== ℓ]
    f_exposures =T_exposures ./ T_transcriptions
    f_exposures_uncoupled = T_exposures_uncoupled ./ T_transcriptions
    push!(fractions_T_exposure, mean(f_exposures))
    push!(std_fractions_T_exposure, std(f_exposures))
    push!(fractions_T_exposure_uncoupled, mean(f_exposures_uncoupled))
    push!(std_fractions_T_exposure_uncoupled, std(f_exposures_uncoupled))
end

plot_df = DataFrame(
    ℓs = ℓs, 
    fractions_T_exposure = -fractions_T_exposure.+1,
    std_fractions_T_exposure = std_fractions_T_exposure,
    mse_fractions_T_exposure = std_fractions_T_exposure/sqrt(1000),
    fractions_T_exposure_uncoupled = -fractions_T_exposure_uncoupled.+1,
    std_fractions_T_exposure_uncoupled = std_fractions_T_exposure_uncoupled,
    mse_fractions_T_exposure_uncoupled = std_fractions_T_exposure_uncoupled/sqrt(1000),
    )
CSV.write("fig/fraction_T_exposure_df_$(label).csv",plot_df)
@pgf t = Table({x = "ℓs", y = "fractions_T_exposure", "y error"="mse_fractions_T_exposure", "col sep"="comma"}, "fraction_T_exposure_df_$(label).csv")
plot_line = Plot(t)
plot_line["mark"]="o"
legend_line = LegendEntry("distance-protected fraction")
push!(axis, plot_line )
push!(axis, legend_line )
@pgf t = Table({x = "ℓs", y = "fractions_T_exposure_uncoupled", "y error"="mse_fractions_T_exposure_uncoupled", "col sep"="comma"}, "fraction_T_exposure_df_$(label).csv")
plot_line = Plot(t)
plot_line["mark"]="*"
legend_line = LegendEntry("coupling-protected fraction")
push!(axis, plot_line )
push!(axis, legend_line )

pgfsave("fig/protected_fraction_plot_$(label).tex",axis)
# pgfsave("fig/exposure_fraction_plot_$(label).svg",axis)


# # plot the mean and std of the mean_distance
# # sort!(df,[:v_translations])
# @pgf axis = Axis(
#     {
#         xlabel = "interaction length \$\\ell\$",
#         ylabel = "distance / codon",
#         grid = "major",
#         "error bars/y dir=both",
#         "error bars/y explicit",
#         legend_pos  = "north west"
#     },
#     # VLine({ dotted, red }, k)
# )
# # group data by rate of translation
# mean_distances = Float64[]
# std_distances = Float64[]

# for ℓ in ℓs
#     distances = df.mean_distance[df.ℓ .== ℓ]
#     push!(mean_distances, mean(distances))
#     push!(std_distances, std(distances))
# end
# plot_line = Plot(Coordinates(ℓs, mean_distances; yerror= std_distances ))
# push!(axis, plot_line)

# pgfsave("fig/mean_distance_plot_$(label).tex",axis)
# pgfsave("fig/mean_distance_plot_$(label).svg",axis)
