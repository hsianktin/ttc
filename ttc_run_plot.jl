using Distributed
using DataFrames
using OffsetArrays
using CSV
using PGFPlotsX
addprocs(7)
label = "ini_pausing"
begin
    using ProgressMeter
end
k = 45
ps = [i for i in 5:70]
L = 1000;
ℓ = 2;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3;
E_c = 1;
k_couple = 100;
k_stalling_0 = 0.3;
k_unstalling_0 = 0.3;
k_ini_pausings = [0.1, 1e2, 1e4];
if isfile("output_$(label).csv")    
    df = CSV.read("output_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        for k_ini_pausing in k_ini_pausings
            cmd = `julia ttc_dynamics.jl $k $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing`
            push!(cmds, cmd)
        end
    end
    @showprogress 1 pmap(run, cmds)

    df = DataFrame(k_couple = Float64[], k_uncouple = Float64[], v_translations = Float64[], v_transcriptions = Float64[], v_stalls = Float64[], v_unstalls = Float64[], k_ini_pausings = Float64[] , L = [], ℓ = [], Eᵦ = Float64[], E_c = Float64[], x₀ = Int[], y₀ = Int[], s₀ = Int[], p₀ = Int[], type=String[], μ = Float64[], σ² = Float64[], μ_c = Float64[])
    # end

    for f in readdir("./data/",join=true)
        temp_df = CSV.read(f,DataFrame)
        for i in 1:length(temp_df[:,1])
            push!(df, temp_df[i,:])
        end
        rm(f)
    end
    CSV.write("output_$(label).csv",df)
end
using LaTeXStrings


# supposing k is the changing variable
# plot df.μ_c as a function of k
sort!(df,[:k_ini_pausings,:v_translations])
@pgf axis = Axis(
    {
        xlabel = L"v_{\mathrm{ribosome}}",
        ylabel = L"p_c",
        grid = "major",
        legend_pos  = "north west"
    },
    VLine({ dotted, red }, k)
)
for k_ini_pausing in k_ini_pausings
    # get the v_translations with the same k_ini_pausing:
    x = df.v_translations[df.k_ini_pausings .== k_ini_pausing]
    y = df.μ_c[df.k_ini_pausings .== k_ini_pausing]
    plot_line = Plot(
        Coordinates(x,y),
    )
    legend_entry = LegendEntry(latexstring("k_{\\mathrm{ini,pause}}=$(k_ini_pausing)"))
    push!(axis,plot_line)
    push!(axis,legend_entry)
end



# to do: multiple plots

pgfsave("fig/line_plot_$(label).tex",axis)
pgfsave("fig/line_plot_$(label).svg",axis)