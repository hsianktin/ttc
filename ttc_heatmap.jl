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
    if isfile("profile/heatmap_$(profile_name).jl")
        include("profile/heatmap_$profile_name.jl")
        println("profile $(profile_name) loaded...")
    else
        error("profile $(profile_name) not found...")
    end
else
    println("Usage: ttc_heatmap.jl [profile_name]")
end

if isfile("output_$(label).csv")
    df = CSV.read("output_$(label).csv",DataFrame)
else
    cmds = []
    for q in qs, p in ps
        cmd = `julia ttc_dynamics.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing`
        push!(cmds, cmd)
    end
    @showprogress 1 pmap(run, cmds)
    # if isfile("output.csv")
    #     df = CSV.read("output.csv",DataFrame)
    # else
    df = DataFrame(
            k_couple = [], 
            k_uncouple = [], 
            v_translations = [], 
            v_transcriptions = [], 
            v_stalls = [], 
            v_unstalls = [], 
            α = Float64[],
            k_ini_pausings = [], 
            L = [], 
            ℓ = [], 
            Eᵦ = [], 
            E_c = [], 
            x₀ = [], 
            y₀ = [], 
            s₀ = [], 
            p₀ = [],
            type=[], 
            μ = [], 
            σ² = [], 
            μ_c = []
        )
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
# using Plots
# backend(:gr)
z = zeros(length(qs),length(ps))
# z = OffsetArray(z,10:90,10:90)

for i in 1:length(df.v_transcriptions)
    x = df.v_translations[i]
    y = df.v_transcriptions[i]
    # if isinteger((x-min)/step+1) && isinteger((y-min)/step+1)
        id_x = round(Int,(x-min)/step+1)
        id_y = round(Int,(y-min)/step+1)
        z[id_x,id_y] = df.μ_c[i]
    # else
    #     println("x = $x, y = $y")
    #     println("(x-min)/step+1 = $((x-min)/step+1), (y-min)/step+1 = $((y-min)/step+1)")    
    # end
end
# using PGFPlotsX to generate heatmap instead of Plots
x = qs
y = ps
@pgf axis = Axis(
    {
        xlabel = "RNAP elongation rate \$q\$",
        ylabel = "ribosome translocation rate \$p\$",
        title = "mean coupling coefficient",
        view = (0, 90),
        colorbar,
        "colormap/jet",
    },
    Plot3(
        {
            surf,
            shader = "flat",
        },
        Coordinates(x, y, transpose(z))
    )
)
# heatmap([i for i in 10:90], [i for i in 10:90],z,size=(500,450))
# xlabel!("v_transcription")
# ylabel!("v_translation")
# title!("mean coupling coefficient coupling_$(Eᵦ)")
plot_df = DataFrame(
    p = [],
    q = [],
    C = []
)
for x in qs, y in ps
    id_x = round(Int,(x-min)/step+1)
    id_y = round(Int,(y-min)/step+1)
        
    push!(plot_df, [
        x,
        y,
        z[id_x,id_y]
    ])
end

CSV.write("fig/mean_coupling_coef_$label.csv",plot_df)
# # pgfsave("fig/mean_coupling_coef_$label.svg",axis)
# pgfsave("fig/mean_coupling_coef_$label.tex",axis)

## plot the heatmap of mean delay/s
z = zeros(Float64,length(qs),length(ps))
for i in 1:length(df.v_transcriptions)
    x = df.v_translations[i]
    y = df.v_transcriptions[i]
    # if isinteger((x-min)/step+1) && isinteger((y-min)/step+1)
        id_x = round(Int,(x-min)/step+1)
        id_y = round(Int,(y-min)/step+1)
        z[id_x,id_y] = df.μ[i]/df.v_translations[i]
    # end
end
# using PGFPlotsX to generate heatmap instead of Plots
x = qs
y = ps
@pgf axis = Axis(
    {
        xlabel = "RNAP elongation rate \$q\$",
        ylabel = "ribosome translocation rate \$p\$",
        title = "mean delay/s",
        view = (0, 90),
        colorbar,
        "colormap/jet",
    },
    Plot3(
        {
            surf,
            shader = "flat",
        },
        Coordinates(x, y, transpose(z))
    )
)

plot_df = DataFrame(
    p = [],
    q = [],
    delay = []
)
for x in qs, y in ps
    id_x = round(Int,(x-min)/step+1)
    id_y = round(Int,(y-min)/step+1)
        
    push!(plot_df, [
        x,
        y,
        z[id_x,id_y]
    ])
end

CSV.write("fig/mean_delay_$label.csv",plot_df)
# # pgfsave("fig/mean_delay_per_second_$label.svg",axis)
# pgfsave("fig/mean_delay_per_second_$label.tex",axis)

## plot the heatmap of std delay/s

z = zeros(Float64,length(qs),length(ps))

for i in 1:length(df.v_transcriptions)
    x = df.v_translations[i]
    y = df.v_transcriptions[i]
    # if isinteger((x-min)/step+1) && isinteger((y-min)/step+1)
        id_x = round(Int,(x-min)/step+1)
        id_y = round(Int,(y-min)/step+1)
        z[id_x,id_y] = sqrt(df.σ²[i])/df.v_translations[i]
    # end
end
# using PGFPlotsX to generate heatmap instead of Plots
x = qs
y = ps
@pgf axis = Axis(
    {
        xlabel = "RNAP elongation rate \$q\$",
        ylabel = "ribosome translocation rate \$p\$",
        title = "standard deviation of delay/s",
        view = (0, 90),
        colorbar,
        "colormap/jet",
    },
    Plot3(
        {
            surf,
            shader = "flat",
        },
        Coordinates(x, y, transpose(z))
    )
)


# pgfsave("fig/std_delay_by_second_$label.svg",axis)
# pgfsave("fig/std_delay_by_second_$label.tex",axis)

# heatmap of protection fraction
if isfile("simu_heatmap_$(label).csv")    
    df = CSV.read("simu_heatmap_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        for q in qs
                cmd = `julia ttc_simu_base.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type`
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
        );    # end

    for f in readdir("./data/",join=true)
        temp_df = CSV.read(f,DataFrame)
        for i in 1:length(temp_df[:,1])
            push!(df, temp_df[i,:])
        end
        # rm(f)
    end
    CSV.write("simu_heatmap_$(label).csv",df)
    for f in readdir("./data/",join=true)
        # temp_df = CSV.read(f,DataFrame)
        # for i in 1:length(temp_df[:,1])
        #     push!(df, temp_df[i,:])
        # end
        rm(f)
    end
end
using LaTeXStrings
using Statistics
z = zeros(Float64,length(qs),length(ps))

for x in ps, y in qs
    temp_df_0 = df[df.v_translations .== x,:]
    temp_df = temp_df_0[temp_df_0.v_transcriptions .== y,:]
    id_x = round(Int,(x-min)/step+1)
    id_y = round(Int,(y-min)/step+1)
    z[id_x,id_y] = mean([1 - temp_df.T_exposure[i]/temp_df.T_transcription[i] for i in 1:length(temp_df.T_transcription)])
end

# using PGFPlotsX to generate heatmap instead of Plots
x = ps
y = qs
@pgf axis = Axis(
    {
        xlabel = "RNAP elongation rate \$q\$",
        ylabel = "ribosome translocation rate \$p\$",
        title = "fraction of protected time \$f_T\$",
        view = (0, 90),
        colorbar,
        "colormap/jet",
    },
    Plot3(
        {
            surf,
            shader = "flat",
        },
        Coordinates(x, y, transpose(z))
    )
)

plot_df = DataFrame(
    q = [],
    p = [],
    f = []
)
for x in qs, y in ps
    id_x = round(Int,(x-min)/step+1)
    id_y = round(Int,(y-min)/step+1)
        
    push!(plot_df, [
        x,
        y,
        z[id_x,id_y]
    ])
end

CSV.write("fig/heatmap_f_T_$label.csv",plot_df)
# pgfsave("fig/heatmap_f_T_$label.tex",axis)
# pgfsave("fig/heatmap_f_T_$label.svg",axis)
df_1 = CSV.read("fig/mean_coupling_coef_$(label).csv",DataFrame)
df_2 = CSV.read("fig/simu_$(label).csv",DataFrame)
df_3 = innerjoin(df_2,df_1, on = [:q,:p])
CSV.write("fig/corr_C_$(label).csv",df_3)
