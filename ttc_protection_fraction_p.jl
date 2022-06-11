using Distributed
using DataFrames
using OffsetArrays
using CSV
using PGFPlotsX
using ProgressMeter
using LaTeXStrings
if length(ARGS) == 1
    profile_name = ARGS[1]
    if isfile("profile/heatmap_$(profile_name).jl")
        include("profile/heatmap_$profile_name.jl")
        println("profile $(profile_name) loaded...")
    else
        error("profile $(profile_name) not found...")
    end
else
    error("Usage: julia ttc_protection_fraction_p.jl [profile_name]")
end

# plot of protection fraction
if isfile("simu_heatmap_$(label).csv")    
    df = CSV.read("simu_heatmap_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        for q in qs
            for i in 1:N
                cmd = `julia ttc_simu_base.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type`
                push!(cmds, cmd)
            end
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
        mean_distance = Float64[]
        );    # end

    for f in readdir("./data/",join=true)
        temp_df = CSV.read(f,DataFrame)
        for i in 1:length(temp_df[:,1])
            push!(df, temp_df[i,:])
        end
        rm(f)
    end
    CSV.write("simu_heatmap_$(label).csv",df)
end
using LaTeXStrings
using Statistics
q = 27
y = q
z = zeros(Float64,length(ps))
δ = zeros(Float64,length(ps))
P = []
F = []
for x in ps
    temp_df_0 = df[df.v_translations .== x,:]
    temp_df = temp_df_0[temp_df_0.v_transcriptions .== y,:]
    id_x = round(Int,(x-min)/step+1)
    z[id_x] = mean([1 - temp_df.T_exposure[i]/temp_df.T_transcription[i] for i in 1:length(temp_df.T_transcription)])
    δ[id_x] = std([1 - temp_df.T_exposure[i]/temp_df.T_transcription[i] for i in 1:length(temp_df.T_transcription)])

end
temp_df = df[df.v_transcriptions .== q,:]
for i in 1:length(temp_df.T_transcription)
    push!(P,temp_df.v_translations[i])
    push!(F,(1-temp_df.T_exposure[i]/temp_df.T_transcription[i]))
end
# using PGFPlotsX to generate heatmap instead of Plots
x = ps
y = qs
@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        xlabel = "ribosome translocation rate \$p\$",
        ylabel = "fraction of protected time \$f_T\$",
        "error bars/y dir=both",
        "error bars/y explicit",
    },
    Plot(
        Coordinates(ps, z; yerror = δ),
    )
)
# save data for plot to CSV file
df_plot = DataFrame(
    p = ps,
    f = z,
    df = δ
)
CSV.write("fig/plot_f_T_p_$(label).csv",df_plot)

pgfsave("fig/plot_f_T_p_$label.tex",axis)
pgfsave("fig/plot_f_T_p_$label.svg",axis)

@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        xlabel = "ribosome translocation rate \$p\$",
        ylabel = "fraction of protected time \$f_T\$",
    },
    Plot(
        {
            only_marks,
            mark_size=1,
            color="blue",
        },
        Coordinates(P,F),
    )
)
# save data for plot to CSV file
df_plot = DataFrame(
    p = P,
    f = F
)
CSV.write("fig/scatter_f_T_p_$(label).csv",df_plot)


pgfsave("fig/scatter_f_T_p_$label.tex",axis)
pgfsave("fig/scatter_f_T_p_$label.svg",axis)

