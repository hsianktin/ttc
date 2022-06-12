mode = "plain"
N = 50000 # repeats
q = 30 # v_transcription
ps = [i for i in 10.0:5:25.0] # v_translation
L = 335;
ℓ = 4;
αs = [10.0^(i) for i in -2:0.25:2];
k_transcription_termination = 0.1;
Eᵦ = 3.0;
E_c = 2.0;
k_couple = 100;
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
if length(ARGS) == 3 
    E_c = parse(Float64, ARGS[3])
end
label = "alphas_var_ps"
type = "kinetic_push"

################
## Simulation ##
################
if isfile("simu_$(label).csv") && !overwrite
    print("simu_$(label).csv exists, skip simulation\n")
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        for α in αs
            cmd = `julia ttc_simu_base.jl $q $p $L $ℓ $α $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $mode $N`
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
sort!(df,[:α,:v_translations])

### Create dataframe 
merged_df = DataFrame(
    α = Float64[],
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
    for α in αs
        T_transcriptions =temp_df.T_transcription[temp_df.α .== α]
        T_translations =temp_df.T_translation[temp_df.α .== α]
        fractions_T_exposure = temp_df.T_exposure[temp_df.α .== α]./temp_df.T_transcription[temp_df.α .== α]
        fractions_T_exposure_uncoupled = temp_df.T_exposure_uncoupled[temp_df.α .== α]./temp_df.T_transcription[temp_df.α .== α]
        push!(merged_df,
            [
                α,
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
push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"\usetikzlibrary{plotmarks}")
markers = ["*","triangle*","square*","pentagon*"]
sizes = [1.8,2,1.8,2]
## effective velocity ##
@pgf axis = SemiLogXAxis(
    {
        xlabel = "translation initiation rate \$\\alpha\$",
        ylabel = "effective transcript. velocity \$\\bar{V}_{\\rm RNAP} \$",
        width = "3.4in",
        height = "3.4in",
    },
)

for i in 1:length(ps)
    t = @pgf Table({x = "α", y = "v_eff_transcriptions", "col sep"="comma"},"merged_df_$(label).csv")
    string = "discard if not={p}{$(ps[i])}"
    t[string]=nothing
    plot_line = @pgf Plot(
        {
            mark=markers[i],
            "mark size"=sizes[i],
        },
        t
    )
    print_tex(t)
    legend_entry = LegendEntry(latexstring("p = $(ps[i])"))
    push!(axis, plot_line)
    push!(axis, legend_entry)
end

pgfsave("fig/effective_velocity_plot_$(label).tex",axis)

## mean fraction of protected time ##
@pgf axis = SemiLogXAxis(
    {
        xlabel = "translation initiation rate \$\\alpha\$",
        ylabel = "mean fraction of protected time",
        width = "3.4in",
        height = "3.4in",
    },
)

for i in 1:length(ps)
    t = @pgf Table({x = "α", y = "fraction_protected", "col sep"="comma"},"merged_df_$(label).csv")
    string = "discard if not={p}{$(ps[i])}"
    t[string]=nothing
    plot_line = @pgf Plot(
        {
            mark=markers[i],
            "mark size"=sizes[i],
        },
        t
    )
    print_tex(t)
    legend_entry = LegendEntry(latexstring("p = $(ps[i])"))
    push!(axis, plot_line)
    push!(axis, legend_entry)
end

pgfsave("fig/fraction_protected_plot_$(label).tex",axis)

## merged ##
@pgf axis = SemiLogXAxis(
    {
        width = "3.4in",
        height = "3.4in",
        ylabel = "ribosome pushing efficiency",
        xlabel = "mean fraction of protected duration",
    },
)
for i in 1:length(ps)
    t = @pgf Table({x = "fraction_protected", y = "efficiency", "col sep"="comma"},"merged_df_$(label).csv")
    string = "discard if not={p}{$(ps[i])}"
    t[string]=nothing
    plot_line = @pgf Plot(
        {
            mark=markers[i],
            "mark size"=sizes[i],
        },
        t
    )
    print_tex(t)
    legend_entry = LegendEntry(latexstring("p = $(ps[i])"))
    push!(axis, plot_line)
    push!(axis, legend_entry)
end

pgfsave("fig/merged_plot_$(label).tex",axis)
