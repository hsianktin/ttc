using DataFrames
using OffsetArrays
using CSV
using PGFPlotsX
using LaTeXStrings

df = DataFrame(k_couple = [], k_uncouple = [], v_translations = [], v_transcriptions = [], v_stalls = [], v_unstalls = [], k_ini_pausings = [], L = [], ℓ = [], Eᵦ = [], E_c = [], x₀ = [], y₀ = [], s₀ = [], p₀ = [],type=[], μ = [], σ² = [], μ_c = [])


if isfile("output_slow_pausing_L_335.csv")
    temp_df = CSV.read("output_slow_pausing_L_335.csv", DataFrame)
    for i in 1:length(temp_df[:,1])
        push!(df, temp_df[i,:])
    end
else
    error("output_slow_pausing_L_335.csv does not exist")
end

# scatter plot instead of heatmap
Z = Float64[]
W = Float64[]
# z = OffsetArray(z,10:90,10:90)

for i in 1:length(df.v_transcriptions)
    x = df.v_translations[i]
    y = df.v_transcriptions[i]
    k_stalling_0 = df.v_stalls[i]
    k_unstalling_0 = df.v_unstalls[i]
    ȳ = y*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)
    push!(W,x/ȳ)
    push!(Z,df.μ_c[i])
end

# save the data of order_parameter=W, coupling_coeff=Z to a csv file in fig/threshold_coupling_coeff.csv
plot_df = DataFrame(order_parameter=W, coupling_coeff=Z)
CSV.write("fig/threshold_coupling_coeff.csv", plot_df)

if isfile("output_low_coupling_L_335.csv")
    temp_df = CSV.read("output_low_coupling_L_335.csv", DataFrame)
    for i in 1:length(temp_df[:,1])
        push!(df, temp_df[i,:])
    end
else
    error("output_low_coupling_L_335.csv does not exist")
end

# scatter plot instead of heatmap
Z = Float64[]
W = Float64[]
# z = OffsetArray(z,10:90,10:90)

for i in 1:length(df.v_transcriptions)
    x = df.v_translations[i]
    y = df.v_transcriptions[i]
    k_stalling_0 = df.v_stalls[i]
    k_unstalling_0 = df.v_unstalls[i]
    ȳ = y*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)
    push!(W,x/ȳ)
    push!(Z,df.μ_c[i])
end

# save the data of order_parameter=W, coupling_coeff=Z to a csv file in fig/threshold_coupling_coeff.csv
plot_df = DataFrame(order_parameter=W, coupling_coeff=Z)
CSV.write("fig/threshold_coupling_coeff_low_coupling.csv", plot_df)



# # using PGFPlotsX to generate heatmap instead of Plots
# t = @pgf Table({x = "order_parameter", y = "coupling_coeff", "col sep"="comma"}, "threshold_coupling_coeff.csv")
# @pgf axis = SemiLogXAxis(
#     {
#         xlabel = "\$p/\\bar q\$",
#         ylabel = "mean coupling coefficient",
#         grid = "major",
#     },
#     Plot(
#         {
#             scatter,
#             "only_marks",
#         },
#         t
#     )
# )
# # heatmap([i for i in 10:90], [i for i in 10:90],z,size=(500,450))
# # xlabel!("v_transcription")
# # ylabel!("v_translation")
# # title!("mean coupling coefficient coupling_$(Eᵦ)")
# pgfsave("fig/threshold_coupling_coeff.tex",axis)

###########################################
##### Fraction of Protecction #############
###########################################

if isfile("simu_heatmap_slow_pausing_L_335.csv")
    df = CSV.read("simu_heatmap_slow_pausing_L_335.csv", DataFrame)
else
    error("simu_heatmap_slow_pausing_L_335.csv does not exist")
end

using LaTeXStrings, Statistics
# scatter plot instead of heatmap
Z = Float64[]
Z₊ = Float64[]  
Z₋ = Float64[]
δZ = Float64[]
W = Float64[]
# z = OffsetArray(z,10:90,10:90)

ps = unique(df.v_translations)
sort!(ps)
qs = unique(df.v_transcriptions)
sort!(qs)
if length(unique(df.v_stalls)) == 1
    k_stalling_0 = unique(df.v_stalls)[1]
else
    error("v_stalls is not unique")
end
if length(unique(df.v_unstalls)) == 1
    k_unstalling_0 = unique(df.v_unstalls)[1]
else
    error("v_unstalls is not unique")
end

for x in [ps[i] for i in 1:length(ps)], y in [qs[i] for i in 1:length(qs)]
    ȳ = y*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)
    temp_df_0 = df[df.v_translations .== x, :]
    temp_df = temp_df_0[temp_df_0.v_transcriptions .== y,:]
    data_array = [1 - temp_df.T_exposure[i]/temp_df.T_transcription[i] for i in 1:length(temp_df.T_transcription)]
    z = mean(data_array)
    δz = std(data_array)
    z₊ = quantile(data_array,0.75) - z
    z₋ = z - quantile(data_array,0.25) 
    push!(W,x/ȳ)
    push!(Z,z)
    push!(Z₊,z₊)
    push!(Z₋,z₋)
    push!(δZ,δz)
end

# save the data of order_parameter=W, coupling_coeff=Z to a csv file in fig/threshold_coupling_coeff.csv
plot_df = DataFrame(order_parameter=W, coupling_coeff=Z, q_plus = Z₊, q_minus = Z₋, δq = δZ)
CSV.write("fig/threshold_f_T.csv", plot_df)


if isfile("simu_heatmap_low_coupling_L_335.csv")
    df = CSV.read("simu_heatmap_low_coupling_L_335.csv", DataFrame)
else
    error("simu_heatmap_low_coupling_L_335.csv does not exist")
end

using LaTeXStrings, Statistics
# scatter plot instead of heatmap
Z = Float64[]
Z₊ = Float64[]  
Z₋ = Float64[]
δZ = Float64[]
W = Float64[]
# z = OffsetArray(z,10:90,10:90)

ps = unique(df.v_translations)
sort!(ps)
qs = unique(df.v_transcriptions)
sort!(qs)
if length(unique(df.v_stalls)) == 1
    k_stalling_0 = unique(df.v_stalls)[1]
else
    error("v_stalls is not unique")
end
if length(unique(df.v_unstalls)) == 1
    k_unstalling_0 = unique(df.v_unstalls)[1]
else
    error("v_unstalls is not unique")
end

for x in [ps[i] for i in 1:length(ps)], y in [qs[i] for i in 1:length(qs)]
    ȳ = y*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)
    temp_df_0 = df[df.v_translations .== x, :]
    temp_df = temp_df_0[temp_df_0.v_transcriptions .== y,:]
    data_array = [1 - temp_df.T_exposure[i]/temp_df.T_transcription[i] for i in 1:length(temp_df.T_transcription)]
    z = mean(data_array)
    δz = std(data_array)
    z₊ = quantile(data_array,0.75) - z
    z₋ = z - quantile(data_array,0.25) 
    push!(W,x/ȳ)
    push!(Z,z)
    push!(Z₊,z₊)
    push!(Z₋,z₋)
    push!(δZ,δz)
end

# save the data of order_parameter=W, coupling_coeff=Z to a csv file in fig/threshold_coupling_coeff.csv
plot_df = DataFrame(order_parameter=W, coupling_coeff=Z, q_plus = Z₊, q_minus = Z₋, δq = δZ)
CSV.write("fig/threshold_f_T_low_coupling.csv", plot_df)


# sort!(plot_df,[:order_parameter])
# t = @pgf Table({x = "order_parameter", y = "coupling_coeff", y_error = "δq", "col sep"="comma"}, "threshold_f_T.csv")
# print_tex(t)
# @pgf axis = SemiLogXAxis(
#     {
#         xlabel = "\$p/\\bar q\$",
#         ylabel = "fraction of protected time \$f_T\$",
#         grid = "major",
#     },
#     Plot(
#         {
#             only_marks,
#             mark = "o",
#             mark_size=1,
#             color="red",
#             # "no marks",
#             "error bars/y dir=both",
#             "error bars/y explicit" ,
#         },
#         t
#     )
# )
# # heatmap([i for i in 10:90], [i for i in 10:90],z,size=(500,450))
# # xlabel!("v_transcription")
# # ylabel!("v_translation")
# # title!("mean coupling coefficient coupling_$(Eᵦ)")
# pgfsave("fig/threshold_f_T.tex",axis)

# ############################
# ### mean vs std ############
# ############################
# t = @pgf Table({x = "coupling_coeff", y = "δq",  "col sep"="comma"}, "threshold_f_T.csv")

# @pgf axis = Axis(
#     {
#         xlabel = "\$\\overline{F}_T\$",
#         ylabel = "\$\\textrm{std} F_T\$",
#         # grid = "major",
#     },
#     Plot(
#         {
#             only_marks,
#             mark = "o",
#             mark_size=1,
#             color="red",
#             # "no marks",
#             # "error bars/y dir=both",
#             # "error bars/y explicit" ,
#         },
#         t
#     )
# )
# pgfsave("fig/mean_std_F_T.tex",axis)

# ############################
# ##### single plot with q = 30
# ############################

# using LaTeXStrings, Statistics
# # scatter plot instead of heatmap
# Z = Float64[]
# Z₊ = Float64[]  
# Z₋ = Float64[]
# W = Float64[]
# # z = OffsetArray(z,10:90,10:90)

# y = q = 30

# for x in ps
#     ȳ = y*(k_unstalling_0)/(k_stalling_0+k_unstalling_0)
#     temp_df_0 = df[df.v_translations .== x, :]
#     temp_df = temp_df_0[temp_df_0.v_transcriptions .== y,:]
#     data_array = [1 - temp_df.T_exposure[i]/temp_df.T_transcription[i] for i in 1:length(temp_df.T_transcription)]
#     z = mean(data_array)
#     z₊ = quantile(data_array,0.75) - z
#     z₋ = z - quantile(data_array,0.25) 
#     push!(W,x)
#     push!(Z,z)
#     push!(Z₊,z₊)
#     push!(Z₋,z₋)
# end
# # save the data of order_parameter=W, coupling_coeff=Z to a csv file in fig/threshold_coupling_coeff.csv
# plot_df = DataFrame(order_parameter=W, coupling_coeff=Z, q_plus = Z₊, q_minus = Z₋)
# CSV.write("fig/scatter_f_T_p.csv", plot_df)

# sort!(plot_df,[:order_parameter])
# t = @pgf Table({x = "order_parameter", y = "coupling_coeff", "col sep"="comma"}, "scatter_f_T_p.csv")
# @pgf axis = Axis(
#     {
#         xlabel = "\$p/\\bar q\$",
#         ylabel = "fraction of protected time \$f_T\$",
#         grid = "major",
#     },
#     Plot(
#         {
#             only_marks,
#             mark="o",
#             mark_size=1,
#             color="red",
#         },
#         Coordinates(plot_df.order_parameter, plot_df.coupling_coeff; yerrorplus=plot_df.q_plus, yerrorminus=plot_df.q_minus)
#     )
# )
# # heatmap([i for i in 10:90], [i for i in 10:90],z,size=(500,450))
# # xlabel!("v_transcription")
# # ylabel!("v_translation")
# # title!("mean coupling coefficient coupling_$(Eᵦ)")
# pgfsave("fig/scatter_f_T_p.tex",axis)
