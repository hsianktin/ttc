mode = "plain"
label = "E_bs_var_ps_long_range" # 
N = 10000 # repeats, deprecated, set to 1 for performance.
q = 30 # v_transcription
ps = [i for i in 5:5:25] # v_translation
L = 335;
ℓ = 40;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦs = [i for i in -2:1:10];
E_c = 2.0;
k_couple = 100;
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
type = "kinetic_push"

################
## Simulation ##
################

if isfile("simu_$(label).csv")   && !overwrite
    print("simu_$(label).csv exists, skip simulation\n")
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        for Eᵦ in Eᵦs
            cmd = `julia ttc_simu_base.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $mode`
            println(cmd)
            push!(cmds, cmd)
        end
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
    # end

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

################
##  Plotting  ##
################


## effective velocity ##
sort!(df,[:Eᵦ,:v_translations])
@pgf axis = Axis(
    {
        xlabel = "binding energy between the ribosome and RNAP \$E_{\\beta}\$",
        ylabel = "effective transcript. velocity \$\\bar{V}_{\\rm RNAP} \$",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north west"
    },
    VLine({ dotted, red }, q*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)),
)

effective_velocity_df = DataFrame(
    Eᵦ = Float64[],
    v_translations = Float64[],
    v_transcriptions = Float64[],
    v_eff_translations = Float64[],
    v_eff_transcriptions = Float64[],
    std_eff_translations = Float64[],
    std_eff_transcriptions = Float64[],
)

for p in ps
    temp_df = df[df.v_translations .== p,:]
    v_eff_transcription = Float64[]
    std_eff_translation = Float64[]
    v_eff_translation = Float64[]
    std_eff_transcription = Float64[]
    for Eᵦ in Eᵦs
        v_eff_transcriptions = (L)./temp_df.T_transcription[temp_df.Eᵦ .== Eᵦ]
        v_eff_translations = (L)./temp_df.T_translation[temp_df.Eᵦ .== Eᵦ]
        push!(v_eff_transcription, mean(v_eff_transcriptions))
        push!(std_eff_translation, std(v_eff_translations))
        push!(v_eff_translation, mean(v_eff_translations))
        push!(std_eff_transcription, std(v_eff_transcriptions))
        push!(effective_velocity_df,
            [
                Eᵦ,
                p,
                q,
                mean(v_eff_translations),
                mean(v_eff_transcriptions),
                std(v_eff_translations),
                std(v_eff_transcriptions)
            ]
            )
    end
    plot_line = Plot(Coordinates(Eᵦs, v_eff_transcription; yerror= std_eff_transcription ))
    legend_line = LegendEntry(latexstring("p=$p"))
    push!(axis, plot_line )
    push!(axis, legend_line )
end

CSV.write(
    "fig/effective_velocity_df_$(label).csv",
    effective_velocity_df
    )

pgfsave("fig/effective_velocity_plot_$(label).tex",axis)
pgfsave("fig/effective_velocity_plot_$(label).svg",axis)
    

## ratio of protection ##
@pgf axis = Axis(
    {
        xlabel = "binding energy between the ribosome and RNAP \$E_{\\beta}\$",
        ylabel = "mean proteccted fraction",
        grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)

fraction_protected_df = DataFrame(
    Eᵦ = Float64[],
    v_translations = Float64[],
    v_transcriptions = Float64[],
    fraction_protected = Float64[],
    std_fraction_exposure = Float64[],
    fraction_protected_uncoupled = Float64[],
    std_fraction_exposure_uncoupled = Float64[],
)

# group data by rate of translation
for p in ps
    temp_df = df[df.v_translations .== p,:]
    fractions_T_exposure = Float64[]
    fractions_T_exposure_uncoupled = Float64[]
    std_fractions_T_exposure = Float64[]
    std_fractions_T_exposure_uncoupled = Float64[]
    for Eᵦ in Eᵦs
        T_exposures = temp_df.T_exposure[temp_df.Eᵦ .== Eᵦ]
        T_transcriptions = temp_df.T_transcription[temp_df.Eᵦ .== Eᵦ]
        T_exposures_uncoupled = temp_df.T_exposure_uncoupled[temp_df.Eᵦ .== Eᵦ]
        f_exposures = T_exposures./T_transcriptions
        f_exposures_uncoupled = T_exposures_uncoupled./T_transcriptions
        push!(fractions_T_exposure, mean(f_exposures))
        push!(std_fractions_T_exposure, std(f_exposures))
        push!(fractions_T_exposure_uncoupled, mean(f_exposures_uncoupled))
        push!(std_fractions_T_exposure_uncoupled, std(f_exposures_uncoupled))
        push!(fraction_protected_df,
            [
                Eᵦ,
                p,
                q,
                1 - mean(f_exposures),
                std(f_exposures),
                1 - mean(f_exposures_uncoupled),
                std(f_exposures_uncoupled)
            ]
        )
    end
    plot_line = Plot(Coordinates(Eᵦs, -fractions_T_exposure .+ 1; yerror= std_fractions_T_exposure ))
    legend_line = LegendEntry(latexstring("p=$p"))
    push!(axis, plot_line )
    push!(axis, legend_line )
end

CSV.write(
    "fig/fraction_protected_df_$(label).csv",
    fraction_protected_df
    )
# savefig in svg and tex format
pgfsave("fig/fraction_protected_plot_$(label).tex",axis)
pgfsave("fig/fraction_protected_plot_$(label).svg",axis)

merged_df = DataFrame(
    Eᵦ = Float64[],
    v_translations = Float64[],
    v_transcriptions = Float64[],
    fraction_protected = Float64[],
    std_fraction_exposure = Float64[],
    v_eff_translations = Float64[],
    v_eff_transcriptions = Float64[],
    std_eff_translations = Float64[],
    std_eff_transcriptions = Float64[],
    efficiency = Float64[],
    std_efficiency = Float64[],
)

for p in ps
    temp_df = df[df.v_translations .== p,:]
    for Eᵦ in Eᵦs
        v_eff_transcriptions = (L)./temp_df.T_transcription[temp_df.Eᵦ .== Eᵦ]
        v_eff_translations = (L)./temp_df.T_translation[temp_df.Eᵦ .== Eᵦ]
        fractions_T_exposure = temp_df.T_exposure[temp_df.Eᵦ .== Eᵦ]./temp_df.T_transcription[temp_df.Eᵦ .== Eᵦ]
        fractions_T_exposure_uncoupled = temp_df.T_exposure_uncoupled[temp_df.Eᵦ .== Eᵦ]./temp_df.T_transcription[temp_df.Eᵦ .== Eᵦ]
        push!(merged_df,
            [
                Eᵦ,
                p,
                q,
                1-mean(fractions_T_exposure),
                std(fractions_T_exposure),
                mean(v_eff_translations),
                mean(v_eff_transcriptions),
                std(v_eff_translations),
                std(v_eff_transcriptions),
                mean(v_eff_translations)/p,
                std(v_eff_translations)/p,
            ]
        )
    end
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

for p in ps
    temp_df = merged_df[merged_df.v_translations .== p,:]
    # plot_line = Plot(Coordinates(temp_df.efficiency, temp_df.fraction_protected; yerror= temp_df.std_fraction_exposure, xerror = temp_df.std_efficiency ))
    plot_line = Plot(Coordinates( temp_df.fraction_protected,temp_df.efficiency))
    legend_line = LegendEntry(latexstring("p=$p"))
    push!(axis, plot_line )
    push!(axis, legend_line )
end
axis
pgfsave("fig/merged_plot_$(label).tex",axis)
pgfsave("fig/merged_plot_$(label).svg",axis)
