mode = "plain"
label = "αs_var_Ls" # 
N = 1000 # repeats, deprecated, set to 1 for performance.
q = 30 # v_transcription
p = 20 # v_translation
Ls = [round(Int, 10.0^k) for k in 2:0.5:6];
ℓ = 4;
αs = [10.0^(k) for k in -3:1:1];
k_transcription_termination = 0.1;
Eᵦ = 3;
E_c = 2;
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
    for α in αs
        for L in Ls
            cmd = `julia ttc_simu_base.jl $q $p $L $ℓ $α $k_transcription_termination $Eᵦ $E_c $k_couple $k_stalling_0 $k_unstalling_0 $k_ini_pausing $d $type $mode $N`
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


include("../ttc_approx.jl")
# group data by rate of translation
alpha = Float64[]
Length = Int[]
v_eff_transcription = Float64[]
std_eff_translation = Float64[]
v_eff_translation = Float64[]
std_eff_transcription = Float64[]
fractions_T_protected = Float64[]
fractions_T_protected_uncoupled = Float64[]
std_fractions_T_protected = Float64[]
std_fractions_T_protected_uncoupled = Float64[]
v_eff_est = Float64[]
F_T_est = Float64[]
C = Float64[]
C₊_est = Float64[]
Cₐ_est = Float64[]
for L in Ls
    temp_df= df[df.L .== L,:]
for α in αs
    S = temp_df.s[temp_df.α .== α]
    C_simu = mean(S)
    T_transcriptions = temp_df.T_transcription[temp_df.α .== α]
    T_translations = temp_df.T_translation[temp_df.α .== α]
    push!(v_eff_transcription, L/mean(T_transcriptions))
    push!(std_eff_transcription, (L/(mean(T_transcriptions)-std(T_transcriptions)) - L/(mean(T_transcriptions)+std(T_transcriptions)) )/2)
    push!(v_eff_translation, L/mean(T_translations))
    push!(std_eff_translation, (L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2)
    T_protecteds = temp_df.T_exposure[temp_df.α .== α]
    T_transcriptions = temp_df.T_transcription[temp_df.α .== α]
    T_protecteds_uncoupled = temp_df.T_exposure_uncoupled[temp_df.α .== α]
    f_protecteds =T_protecteds ./ T_transcriptions
    f_protecteds_uncoupled = T_protecteds_uncoupled ./ T_transcriptions
    push!(fractions_T_protected, 1-mean(f_protecteds))
    push!(std_fractions_T_protected, std(f_protecteds))
    push!(fractions_T_protected_uncoupled, mean(f_protecteds_uncoupled))
    push!(std_fractions_T_protected_uncoupled, std(f_protecteds_uncoupled))
    push!(v_eff_est, V̄ₐ(p,q,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c,ℓ))
    push!(F_T_est, Fₜ₊(p,q,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c,ℓ,27))
    push!(C₊_est, C₊(p,q,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c))
    push!(Cₐ_est, C₊(p,q,k_couple,k_unstalling_0,k_stalling_0,Eᵦ,E_c))
    push!(C, C_simu)
    push!(alpha,α)
    push!(Length,L)
end
end
plot_df = DataFrame(
        α=alpha,
        L=Length, 
        v_eff_transcription=v_eff_transcription, 
        std_eff_transcription=std_eff_transcription, 
        v_eff_translation=v_eff_translation, 
        std_eff_translation=std_eff_translation, 
        fractions_T_protected=fractions_T_protected,
        std_fractions_T_protected=std_fractions_T_protected,
        fractions_T_protected_uncoupled=fractions_T_protected_uncoupled,
        std_fractions_T_protected_uncoupled=std_fractions_T_protected_uncoupled,
        v_eff_est=v_eff_est,
        F_T_est=F_T_est,
        C₊_est = C₊_est,
        Cₐ_est = Cₐ_est,
        C = C,
    )
CSV.write("fig/simu_df_$(label).csv",plot_df)







################
##  Plotting  ##
################


# ## effective velocity ##
# push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"\pgfplotsset{
#     discard if not/.style 2 args={
#         x filter/.code={
#             \edef\tempa{\thisrow{#1}}
#             \edef\tempb{#2}
#             \ifx\tempa\tempb
#             \else
#                 \def\pgfmathresult{inf}
#             \fi
#         }
#     }
# }")
# #### Creating DataFrame ####

# merged_df = DataFrame(
#     L = Float64[],
#     α = Float64[],
#     v_transcriptions = Float64[],
#     fraction_protected = Float64[],
#     std_fraction_exposure = Float64[],
#     v_eff_translations = Float64[],
#     v_eff_transcriptions = Float64[],
#     std_eff_translations = Float64[],
#     std_eff_transcriptions = Float64[],
#     efficiency = Float64[],
#     std_efficiency = Float64[],
# )

# for L in Ls
#     temp_df = df[df.L .== L,:]
#     for α in αs
#         T_transcriptions = temp_df.T_transcription[temp_df.α .== α]
#         T_translations = temp_df.T_translation[temp_df.α .== α]
#         fractions_T_exposure = temp_df.T_exposure[temp_df.α .== α]./temp_df.T_transcription[temp_df.α .== α]
#         fractions_T_exposure_uncoupled = temp_df.T_exposure_uncoupled[temp_df.α .== α]./temp_df.T_transcription[temp_df.α .== α]
#         push!(merged_df,
#             [
#                 L,
#                 α,
#                 q,
#                 1-mean(fractions_T_exposure),
#                 std(fractions_T_exposure),
#                 L/mean(T_translations),
#                 L/mean(T_transcriptions),
#                 (L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2,
#                 (L/(mean(T_transcriptions)-std(T_transcriptions)) - L/(mean(T_transcriptions)+std(T_transcriptions)) )/2,
#                 (L/mean(T_translations))/α,
#                 ((L/(mean(T_translations)-std(T_translations)) - L/(mean(T_translations)+std(T_translations)) )/2)/p,
#             ]
#         )
#     end
# end

# CSV.write(
#     "fig/merged_df_$(label).csv",
#     merged_df
# )

# ## effective velocity ##
# sort!(df,[:L,:v_translations])
# @pgf axis = Axis(
#     {
#         width = "3in",
#         height = "3in",
#         clip = "false",
#         xlabel = "ribosome translocation rate \$α\$",
#         ylabel = "effective transcript. velocity \$\\bar{V}_{\\rm RNAP} \$",
#         # grid = "major",
#         "error bars/y dir=both",
#         "error bars/y explicit",
#         legend_pos  = "north east"
#     },
#     # VLine({ dotted, red }, q*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)),
# )


# for L in Ls
#     t = @pgf Table({x="α",y="v_eff_transcriptions", "col sep"="comma"},"merged_df_$(label).csv")
#     t["discard if not={L}{$(L)}"]=nothing
#     plot_line = Plot(t)
#     legend_line = LegendEntry(latexstring("L=$(L)"))
#     push!(axis, plot_line )
#     push!(axis, legend_line )
# end

# pgfsave("fig/effective_velocity_plot_$(label).tex",axis)
# # pgfsave("fig/effective_velocity_plot_$(label).svg",axis)
    

# ## ratio of protection ##
# @pgf axis = Axis(
#     {        
#         width = "3in",
#         height = "3in",
#         clip = "false",
#         xlabel = "ribosome translocation rate \$α\$",
#         ylabel = "mean fraction of protection \$ F_T\$",
#         legend_pos  = "north west"
#     },
#     # VLine({ dotted, red }, q)
# )

# for L in Ls
#     t = @pgf Table({x="α",y="fraction_protected", "col sep"="comma"},"merged_df_$(label).csv")
#     t["discard if not={L}{$(L)}"]=nothing
#     plot_line = Plot(t)
#     legend_line = LegendEntry(latexstring("L=$(L)"))
#     push!(axis, plot_line )
#     push!(axis, legend_line )
# end
# # savefig in svg and tex format
# pgfsave("fig/fraction_exposure_plot_$(label).tex",axis)
# # pgfsave("fig/fraction_exposure_plot_$(label).svg",axis)

# @pgf axis = Axis(
#     {
#         width = "3in",
#         height = "3in",
#         clip = "false",
#         ylabel = "ribosome pushing efficiency \$ \\eta = \\bar V_{\\rm rib}/α\$",
#         xlabel = "mean protected fraction \$ F_T\$",
#         colorbar,
#         "colorbar style"=@pgf {width="0.2cm"}
#         # "error bars/y dir=both",
#         # "error bars/y explicit",
#         # "error bars/x dir=both",
#         # "error bars/x explicit",
#     },
# )

# # for L in Ls
# #     temp_df = merged_df[merged_df.L .== L,:]
# #     # plot_line = Plot(Coordinates(temp_df.efficiency, temp_df.fraction_protected; yerror= temp_df.std_fraction_exposure, xerror = temp_df.std_efficiency ))
# #     plot_line = Plot(Coordinates( temp_df.fraction_protected,temp_df.efficiency))
# #     legend_line = LegendEntry(latexstring("L=$(L)"))
# #     push!(axis, plot_line )
# #     push!(axis, legend_line )
# # end

# for L in Ls
#     t = @pgf Table({x="fraction_protected",y="efficiency", "col sep"="comma", "meta"="α"},"merged_df_$(label).csv")
#     t["discard if not={L}{$(L)}"]=nothing
#     plot_line = @pgf PlotInc({"scatter", "scatter src"="explicit"},t)
#     legend_line = LegendEntry(latexstring("L=$(L)"))
#     push!(axis, plot_line )
#     push!(axis, legend_line )
# end


# # axis
# pgfsave("fig/merged_plot_$(label).tex",axis)
# # pgfsave("fig/merged_plot_$(label).svg",axis)