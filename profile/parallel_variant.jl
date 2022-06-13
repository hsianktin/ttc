modes = ["plain","profile₁","profile₂"]
N = 10000 # 10s of repeats. The total number is 100.
q = 30 # v_transcription
ps = [i for i in 3:2:30] # v_translation
L = 335;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3.0;
E_c = 2.0;
k_couple = 100;
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
type = "kinetic_push"
if length(ARGS) == 3 
    E_c = parse(Float64, ARGS[3])
end
label = "mode_variant_E_c_$(E_c)_E_b_$(Eᵦ)_alpha_$(k_translation_initiation)" # 


if isfile("simu_$(label).csv")   && !overwrite
    print("simu_$(label).csv exists, skip simulation\n")
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        for mode in modes
            cmd = `julia ttc_simu_base.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $mode $N`
            # println(cmd)
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
merged_df = DataFrame(
    mode = String[],
    p = Float64[],
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

for mode in modes
    temp_df = df[df.mode .== mode,:]
    for p in ps
        T_transcriptions = temp_df.T_transcription[temp_df.v_translations .== p]
        T_translations = temp_df.T_translation[temp_df.v_translations .== p]
        fraction_T_protected = -temp_df.T_exposure[temp_df.v_translations .== p]./temp_df.T_transcription[temp_df.v_translations .== p] .+1
        fraction_T_protected_uncoupled = - temp_df.T_exposure_uncoupled[temp_df.v_translations .== p]./temp_df.T_transcription[temp_df.v_translations .== p] .+1
        push!(merged_df,
            [
                mode,
                p,
                q,
                mean(fraction_T_protected),
                std(fraction_T_protected),
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
##### Specify Plotting Arguments
using LaTeXStrings
sort!(df,[:mode,:v_translations])
@pgf axis = Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "ribosome translocation rate \$p\$",
        ylabel = "effective transcript. velocity \$\\bar{V}_{\\rm RNAP} \$",
        # grid = "major",
        # "error bars/y dir=both",
        # "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, q*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)),
)

function tex_mode(mode)
    if mode == "plain"
        return "uniform"
    elseif mode == "profile₁"
        return "profile 1"
    elseif mode == "profile₂"
        return "profile 2"
    end
end

for mode in modes
    t = @pgf Table({x = "p", y="v_eff_transcriptions", "col sep"="comma"}, "merged_df_$(label).csv")
    t["discard if not={mode}{$mode}"]=nothing
    print_tex(t)
    println()
    plot_line = Plot(t)
    legend_line = LegendEntry(tex_mode(mode))
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
        ylabel = "mean protected fraction \$F_T\$",
        # grid = "major",
        # "error bars/y dir=both",
        # "error bars/y explicit",
        legend_pos  = "north west"
    },
    # VLine({ dotted, red }, k)
)

for mode in modes
    t = @pgf Table({x = "v_translations", y="fraction_protected", "col sep"="comma"}, "merged_df_$(label).csv")
    t["discard if not={mode}{$mode}"]=nothing
    print_tex(t)
    println()
    plot_line = Plot(t)
    legend_line = LegendEntry(tex_mode(mode))
    push!(axis, plot_line )
    push!(axis, legend_line )
end

# savefig in tex format
pgfsave("fig/fraction_protected_plot_$(label).tex",axis)


#### merged_plot
@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        ylabel = "ribosome pushing efficiency \$ \\eta = \\bar V_{\\rm rib}/p\$",
        xlabel = "mean protected fraction \$F_T\$",
        # "error bars/y dir=both",
        # "error bars/y explicit",
        # "error bars/x dir=both",
        # "error bars/x explicit",
    },
)

for mode in modes
    t = @pgf Table({x = "fraction_protected", y="efficiency", "col sep"="comma"}, "merged_df_$(label).csv")
    t["discard if not={mode}{$mode}"]=nothing
    print_tex(t)
    println()
    plot_line = Plot(t)
    legend_line = LegendEntry(tex_mode(mode))
    push!(axis, plot_line )
    push!(axis, legend_line )
end

pgfsave("fig/merged_plot_$(label).tex",axis)
