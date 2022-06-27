using Distributed
using DataFrames
using OffsetArrays
using CSV
using PGFPlotsX
using LaTeXStrings
addprocs(8)
begin
    using ProgressMeter
end

if length(ARGS) == 1
    profile_name = ARGS[1]
    overwrite = false
    
elseif length(ARGS) == 2
    profile_name = ARGS[1]
    overwrite = parse(Bool,ARGS[2])
else
    error("Usage: ttc_delay_bifurcation.jl [profile_name]")
end
if isfile("profile/delay_bifurcation_$(profile_name).jl")
    include("profile/delay_bifurcation_$profile_name.jl")
    println("profile $(profile_name) loaded...")
else
    error("profile $(profile_name) not found...")
end
# just reusing the code from ttc_heatmap.jl, but don't really need the `df`
if isfile("fig/delay_$(label).csv") && !overwrite
    plot_df = CSV.read("fig/delay_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        cmd = `julia ttc_dynamics.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $N`
        push!(cmds, cmd)
    end
    @showprogress 1 pmap(run, cmds)
    df = DataFrame(k_couple = [], k_uncouple = [], v_translations = [], v_transcriptions = [], v_stalls = [], v_unstalls = [], k_ini_pausings = [], L = [], ℓ = [], Eᵦ = [], E_c = [], x₀ = [], y₀ = [], s₀ = [], p₀ = [],type=[], μ = [], σ² = [], μ_c = [])
    for f in readdir("./data/",join=true)
        temp_df = CSV.read(f,DataFrame)
        for i in 1:length(temp_df[:,1])
            push!(df, temp_df[i,:])
        end
        rm(f)
    end
end

t_linspace = [i for i in .1:.1:20]
Y = length(t_linspace)
X = length(ps)
z = zeros(X,Y)
plot_df = DataFrame(t = [], p = [], pdf = [], pdf_relat = [])
for x in 1:X
    v_transcription_0 = q
    v_translation_0 = ps[x]
    temp_df = CSV.read("fig/delay/pdf_$(v_transcription_0)_$(v_translation_0)_$(k_stalling_0)_$(k_unstalling_0)_$(k_couple)_$(k_uncouple)_$(E_c)_$(Eᵦ).csv",DataFrame)
    for y in 1:Y
        z[x,y] = (temp_df.pdf_linspace[y])/maximum(temp_df.pdf_linspace)
        push!(plot_df, [t_linspace[y], ps[x],(temp_df.pdf_linspace[y]) ,z[x,y]])
    end
end

CSV.write("fig/delay_$(label).csv",plot_df)
t = @pgf Table({
    x = "t",
    y = "p",
    z = "pdf_relat",
    "col sep" = "comma",
},"delay_$(label).csv")
print_tex(t)
@pgf axis = Axis(
    {
        width = "3.4in",
        height = "3.4in",
        ylabel = "RNAP elongation rate \$q\$",
        xlabel = "delay time / s",
        title = "delay time distribution",
        view = (0, 90),
        colorbar,
        "colormap/hot2",
        "point meta max=1, point meta min=0",
    },
    Plot3(
        {
            surf,
            shader = "flat",
            "mesh/rows" = X,
            "mesh/cols" = Y,
        },
        t
    )
)

pgfsave("fig/delay_$(label).tex",axis)
# pgfsave("fig/delay_$(label).svg",axis)
