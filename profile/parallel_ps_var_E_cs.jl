mode = "plain"
label = "ps_var_E_cs" # 
N = 50000 # repeats, deprecated, set to 1 for performance.
q = 30 # v_transcription
ps = [i for i in 2:2:60] # v_translation
L = 335;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3;
E_cs = [i for i in 0:0.5:2];
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
        for E_c in E_cs
            cmd = `julia ttc_simu_base.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $mode $N`
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

merged_df = DataFrame(
    Ec = Float64[],
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

for E_c in E_cs
    temp_df = df[df.E_c .== E_c,:]
    for p in ps
        T_transcriptions = temp_df.T_transcription[temp_df.v_translations .== p]
        T_translations = temp_df.T_translation[temp_df.v_translations .== p]
        fractions_T_exposure = temp_df.T_exposure[temp_df.v_translations .== p]./temp_df.T_transcription[temp_df.v_translations .== p]
        fractions_T_exposure_uncoupled = temp_df.T_exposure_uncoupled[temp_df.v_translations .== p]./temp_df.T_transcription[temp_df.v_translations .== p]
        push!(merged_df,
            [
                E_c,
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
sort!(df,[:E_c,:v_translations])
@pgf axis = Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "ribosome translocation rate \$p\$",
        ylabel = "effective transcript. velocity \$\\bar{V}_{\\rm RNAP} \$",
        # grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "north east"
    },
    # VLine({ dotted, red }, q*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)),
)


for E_c in E_cs
    t = @pgf Table({x="p",y="v_eff_transcriptions", "col sep"="comma"},"merged_df_$(label).csv")
    t["discard if not={Ec}{$(E_c)}"]=nothing
    plot_line = Plot(t)
    legend_line = LegendEntry(latexstring("E_+=$(E_c)"))
    push!(axis, plot_line )
    push!(axis, legend_line )
end

pgfsave("fig/effective_velocity_plot_$(label).tex",axis)
# pgfsave("fig/effective_velocity_plot_$(label).svg",axis)
    

## ratio of protection ##
@pgf axis = Axis(
    {        
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "ribosome translocation rate \$p\$",
        ylabel = "mean fraction of protection \$ F_T\$",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)

for E_c in E_cs
    t = @pgf Table({x="p",y="fraction_protected", "col sep"="comma"},"merged_df_$(label).csv")
    t["discard if not={Ec}{$(E_c)}"]=nothing
    plot_line = Plot(t)
    legend_line = LegendEntry(latexstring("E_+=$(E_c)"))
    push!(axis, plot_line )
    push!(axis, legend_line )
end
# savefig in svg and tex format
pgfsave("fig/fraction_exposure_plot_$(label).tex",axis)
# pgfsave("fig/fraction_exposure_plot_$(label).svg",axis)

@pgf axis = Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        ylabel = "ribosome pushing efficiency \$ \\eta = \\bar V_{\\rm rib}/p\$",
        xlabel = "mean protected fraction \$ F_T\$",
        colorbar,
        "colorbar style"=@pgf {width="0.2cm"}
        # "error bars/y dir=both",
        # "error bars/y explicit",
        # "error bars/x dir=both",
        # "error bars/x explicit",
    },
)

# for E_c in E_cs
#     temp_df = merged_df[merged_df.E_c .== E_c,:]
#     # plot_line = Plot(Coordinates(temp_df.efficiency, temp_df.fraction_protected; yerror= temp_df.std_fraction_exposure, xerror = temp_df.std_efficiency ))
#     plot_line = Plot(Coordinates( temp_df.fraction_protected,temp_df.efficiency))
#     legend_line = LegendEntry(latexstring("E_c=$(E_c)"))
#     push!(axis, plot_line )
#     push!(axis, legend_line )
# end

for E_c in E_cs
    t = @pgf Table({x="fraction_protected",y="efficiency", "col sep"="comma", "meta"="p"},"merged_df_$(label).csv")
    t["discard if not={Ec}{$(E_c)}"]=nothing
    plot_line = @pgf PlotInc({"scatter", "scatter src"="explicit"},t)
    legend_line = LegendEntry(latexstring("E_+=$(E_c)"))
    push!(axis, plot_line )
    push!(axis, legend_line )
end


# axis
pgfsave("fig/merged_plot_$(label).tex",axis)
# pgfsave("fig/merged_plot_$(label).svg",axis)