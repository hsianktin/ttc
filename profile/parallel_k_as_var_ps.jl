mode = "plain"
label = "k_as_var_ps" # 
N = 1000 # repeats, deprecated, set to 1 for performance.
q = 30 # v_transcription
ps = [i for i in 10:5.0:25] # v_translation
L = 335;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3.0;
E_c = 2.0;
k_couples = [10.0^i for i in 0:0.2:3];
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
type = "kinetic_push"
# overwrite = true

################
## Simulation ##
################

if isfile("simu_$(label).csv") && !overwrite
    print("simu_$(label).csv exists, skip simulation\n")
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        for k_couple in k_couples
            cmd = `julia ttc_simu_base.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $mode $N`
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

################
##  Plotting  ##
################
push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"\pgfplotsset{
    discard if not/.style 2 args={
        x filter/.code={
            \edef\tempa{\thisrow{#1}}
            \edef\tempb{#2}
            \ifx\tempa\tempb
            \else
                \def\pgfmathresult{inf}
            \fi
        }
    }
}")
#### Creating DataFrame ####
sort!(df,[:k_couple,:v_translations])
merged_df = DataFrame(
    k_couple = Float64[],
    p = Float64[],
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
    for k_couple in k_couples
        T_transcriptions = temp_df.T_transcription[temp_df.k_couple .== k_couple]
        T_translations = temp_df.T_translation[temp_df.k_couple .== k_couple]
        fractions_T_exposure = temp_df.T_exposure[temp_df.k_couple .== k_couple]./temp_df.T_transcription[temp_df.k_couple .== k_couple]
        fractions_T_exposure_uncoupled = temp_df.T_exposure_uncoupled[temp_df.k_couple .== k_couple]./temp_df.T_transcription[temp_df.k_couple .== k_couple]
        push!(merged_df,
            [
                k_couple,
                p,
                q,
                1-mean(fractions_T_exposure),
                std(fractions_T_exposure),
                L/mean(T_translations),
                L/mean(T_transcriptions),
                (L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2,
                (L/(mean(T_transcriptions)-std(T_transcriptions)) - L/(mean(T_transcriptions)+std(T_transcriptions)) )/2,
                (L/mean(T_translations))/p,
                ((L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2)/p,
            ]
        )
    end
end


CSV.write(
    "fig/merged_df_$(label).csv",
    merged_df
)


## effective velocity ##
@pgf axis = SemiLogXAxis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "ribosome-RNAP coupling rate \$k_{\\textrm{a}}\$",
        ylabel = "effective transcript. velocity \$\\bar{V}_{\\rm RNAP} \$",
        # grid = "major",
        # "error bars/y dir=both",
        # "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, q*(v_unstalling_0)/(v_stalling_0+v_unstalling_0)),
)

for p in ps
    t = @pgf Table({x = "k_couple", y = "v_eff_transcriptions", "col sep"="comma"}, "merged_df_$(label).csv")
    t["discard if not={p}{$(p)}"]=nothing
    print_tex(t)
    println()
    plot_line = Plot(t)
    legend_line = LegendEntry("\$p = $(p)\$")
    push!(axis, plot_line )
    push!(axis, legend_line )
end
pgfsave("fig/$(label)_effective_velocity_plot.tex",axis)
    

## ratio of protection ##
@pgf axis = SemiLogXAxis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "ribosome-RNAP coupling rate \$k_{\\textrm{a}}\$",
        ylabel = "mean proteccted fraction \$F_T \$",
        # grid = "major",
        # "error bars/y dir=both",
        # "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)
for p in ps
    t = @pgf Table({x = "k_couple", y = "fraction_protected", "col sep"="comma"}, "merged_df_$(label).csv")
    t["discard if not={p}{$(p)}"]=nothing
    print_tex(t)
    println()
    plot_line = Plot(t)
    legend_line = LegendEntry("\$p = $(p)\$")
    push!(axis, plot_line )
    push!(axis, legend_line )
end

# savefig in tex format
pgfsave("fig/$(label)_fraction_protected_plot.tex",axis)

@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        clip = "false",
        ylabel = "ribosome pushing efficiency \$ \\eta = \\bar V_{\\rm rib}/p\$",
        xlabel = "mean proteccted fraction \$F_T \$",
        # "error bars/y dir=both",
        # "error bars/y explicit",
        # "error bars/x dir=both",
        # "error bars/x explicit",
        legend_pos  = "north west"
    },
)

for p in ps
    t = @pgf Table({x = "fraction_protected", y = "efficiency", "col sep"="comma"}, "merged_df_$(label).csv")
    t["discard if not={p}{$(p)}"]=nothing
    print_tex(t)
    println()
    plot_line = Plot(t)
    legend_line = LegendEntry("\$p = $(p)\$")
    push!(axis, plot_line )
    push!(axis, legend_line )
end
# axis
pgfsave("fig/$(label)_merged_plot.tex",axis)
