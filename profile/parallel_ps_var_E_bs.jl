profile = "plain"
label = "ps_var_E_bs" # 
N = 10000 # repeats, deprecated, set to 1 for performance.
q = 30 # v_transcription
ps = [i for i in 5:1:25] # v_translation
L = 335;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦs = [-2.5,0,2.5,5,7.5];
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

if isfile("simu_$(label).csv") && !overwrite
    print("simu_$(label).csv exists, skip simulation\n")
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        for Eᵦ in Eᵦs
            cmd = `julia ttc_simu_base.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $profile $N`
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

using StatsBase, LinearAlgebra
using Pipe

function pdf_corr(X,Y)
    ∑(x) = sum(x)
    min = @pipe minimum([X;Y]) |> floor(Int,_*10)/10
    max = @pipe maximum([X;Y]) |> ceil(Int,_*10)/10
    dx = (max - min)/minimum([length(X),length(Y)])*10
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
df_corr = @pipe df |> 
groupby(_, [:v_translations, :Eᵦ]) |>
combine(_,
    [:T_transcription, :T_translation] => (
        (X,Y) -> (corr = pdf_corr(X,Y))
        ) => :corr 
    ) |>  
    sort(_, [:v_translations, :Eᵦ]) |> rename(_, :v_translations => :p)
CSV.write(
    "fig/pdf_corr_df_$(label).csv",
    df_corr
    )
@pgf axis = Axis(
    {
        width = "3in",
        height = "3in",
        clip = "false",
        xlabel = "ribosome velocity \$p\$",
        ylabel = "correlation",
        # grid = "major",
        "error bars/y dir=both",
        "error bars/y explicit",
        legend_pos  = "south east"
    },
    # VLine({ dotted, red }, q*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)),
)

for Eᵦ in Eᵦs
    t = @pgf Table({x = "p", y = "corr", "col sep"="comma"}, "pdf_corr_df_$(label).csv")
    t["discard if not={Eᵦ}{$(convert(Float64,Eᵦ))}"]=nothing
    plot_line = Plot(t)
    legend_line = LegendEntry("\$E_{a}\$=$Eᵦ")
    push!(axis, plot_line)
    push!(axis, legend_line)
end

pgfsave("fig/pdf_corr_plot_$(label).tex",axis)
        