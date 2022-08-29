profile = "plain"
min = 3
step = 0.5
max = 30
qs = [k for k in min:step:max]
ps = [p for p in min:step:max]
L = 335;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3;
E_c = 2;
k_couple = 100;
k_stalling_0 = 3;
k_unstalling_0 = 3;
k_ini_pausing = k_stalling_0;
d = 27;
N = 1000
label = "fast_pausing_L_$(L)"
type = "kinetic_push"

################
## Simulation ##
################
overwrite = false
if isfile("simu_$(label).csv") && !overwrite
    print("simu_$(label).csv exists, skip simulation\n")
    df = CSV.read("simu_$(label).csv",DataFrame)
else
    cmds = []
    for p in ps
        for q in qs
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

using StatsBase, LinearAlgebra
using Pipe

function pdf_corr(X,Y)
    ∑(x) = sum(x)
    min = @pipe minimum([X;Y]) |> floor(Int,_*10)/10
    max = @pipe maximum([X;Y]) |> ceil(Int,_*10)/10
    dx = (max - min)/minimum([length(X),length(Y)])*√(minimum([length(X),length(Y)]))
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

# group df by the value of ps and qs
df_result = @pipe df |> 
    groupby(_, [:v_translations, :v_transcriptions]) |>
    combine(_,
        [:T_transcription, :T_translation] => (
            (X,Y) -> (corr = pdf_corr(X,Y))
            ) => :corr 
    ) |>
    sort(_, [:v_translations, :v_transcriptions]) |> 
    rename(_, :v_translations => :p) |> rename(_, :v_transcriptions => :q)


using CSV
CSV.write("./fig/simu_$(label).csv",df_result)

# using PGFPlotsX
# t = @pgf Table({x = "v_translations", y = "v_transcriptions", z = "corr", "col sep" = "comma"}, "simu_$(label).csv")
# @pgf axis = Axis(
#     {
#         width = "3.4in",
#         height = "3.4in",
#         xlabel = "death rate \$\\mu\$",
#         ylabel = "birth rate \$\\beta\$",
#         title = "population at equilibrium \$ N \$",
#         legend_pos="south east",
#         view = (0, 90),
#         colorbar,
#         "colorbar style"=@pgf {width="0.2cm"}
#     },
#     Plot3(
#         {
#             surf,
#             shader = "flat",
#             "mesh/rows" = length(unique(df_result.v_transcriptions)),
#             "mesh/cols" = length(unique(df_result.v_translations)),
#         },
#         t
#     ),
#     # LegendEntry(type2latex(type))
# )
# pgfsave("./fig/simu_$(label).tex", axis)