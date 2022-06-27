# for a simple model, we can calculate the probability density explicitly.
using DataFrames
using CSV
using OffsetArrays
using Statistics
using PGFPlotsX
using LaTeXStrings
## initialize parameters
v_transcription_0 = 70.0
v_translation_0 = 40.0
L = 1000 # length of the sequence
ℓ = 2 # characteristic range of couping interactions
k_translation_initiation = 1.0 # translation initiation rate
k_transcription_termination = 0.1 # transcription termination rate
Eᵦ = 3 # energy cost of uncoupling
E_c = 1 # energy contributing to pushing RNAP out of stalling state
# v_y′ = v_y * exp(-Eᵦ)
x₀ = 0 # initial position of the ribosome
y₀ = 1 # initial position of the RNAP
s₀ = 0 # initially uncoupled state
p₀ = 0 # initially processive RNAP
k_couple = 100 # rate of coupling
k_uncouple = 100*exp(-Eᵦ) # rate of uncoupling
type = "kinetic_push"
# velocity_landscape = DataFrame(translation = Float32[],transcription = Float32[]);
# stalling_landscape = DataFrame(stalling = Float32[],unstalling = Float32[]);

k_stalling_0 = 0.2 # initial stalling rate
k_unstalling_0 = 0.2 # initial unstalling rate
k_ini_pausing = 1e3
# accept ARGS as parameters
if length(ARGS) == 12
    v_transcription_0 = parse(Float64,ARGS[1])
    v_translation_0 = parse(Float64,ARGS[2])
    L = parse(Int64,ARGS[3])
    ℓ = parse(Int64,ARGS[4])
    k_translation_initiation = parse(Float64,ARGS[5])
    k_transcription_termination = parse(Float64,ARGS[6])
    Eᵦ = parse(Float64,ARGS[7])
    E_c = parse(Float64,ARGS[8])
    k_couple = parse(Float64,ARGS[9])
    k_uncouple = k_couple * exp(-Eᵦ)
    k_stalling_0 = parse(Float64,ARGS[10])
    k_unstalling_0 = parse(Float64,ARGS[11])
    k_ini_pausing = parse(Float64,ARGS[12])
    plot_flag = 0
elseif length(ARGS) == 13
    v_transcription_0 = parse(Float64,ARGS[1])
    v_translation_0 = parse(Float64,ARGS[2])
    L = parse(Int64,ARGS[3])
    ℓ = parse(Int64,ARGS[4])
    k_translation_initiation = parse(Float64,ARGS[5])
    k_transcription_termination = parse(Float64,ARGS[6])
    Eᵦ = parse(Float64,ARGS[7])
    E_c = parse(Float64,ARGS[8])
    k_couple = parse(Float64,ARGS[9])
    k_uncouple = k_couple * exp(-Eᵦ)
    k_stalling_0 = parse(Float64,ARGS[10])
    k_unstalling_0 = parse(Float64,ARGS[11])
    k_ini_pausing = parse(Float64,ARGS[12])
    plot_flag = parse(Int,ARGS[13])
else
    print("incorrect number of arguments, use default values")
end

# println("velocity_landscape.csv not found, initialize a uniform landscape based on base rates")
velocity_landscape = DataFrame(translation = [v_translation_0 for i in 0:L],transcription = [v_transcription_0 for i in 1:L+1]);
velocity_landscape.translation[1] = k_translation_initiation # the first step is nothing but the initiation step
velocity_landscape.transcription[L+1] = k_transcription_termination # the last step is nothing but the termination step
# println("stalling_landscape.csv not found")
stalling_landscape = DataFrame(stalling = [k_stalling_0 for i in 1:L+1],unstalling = [k_unstalling_0 for i in 1:L+1]);
# "deterministic" stalling at certain position [[@Hatoum2008]] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2819085/
stalling_landscape.stalling[y₀] = k_ini_pausing
stalling_landscape.unstalling[y₀] = 0.01
v_translations = velocity_landscape.translation # translation rate nt/second
v_translations = OffsetArray(v_translations,0:L) # reset the array, such that v_translations[x] represents the translation elongation or initiation rate at position x
v_transcriptions = velocity_landscape.transcription # transcription rate nt/second
v_stalls = stalling_landscape.stalling # stall rate
v_unstalls = stalling_landscape.unstalling # unstall rate

## initialize probability distribution assuming known initial x₀ and y₀
P₀ = zeros(Float64,y₀+1,2,2)
P₀ = OffsetArray(P₀,0:y₀,0:1,0:1) # the 0 in the first dimension represents the absence of RNAP
P₀[x₀,s₀,p₀] = 1
P = P₀

include("ttc_utils.jl")

global P = P₀
for i = y₀ : L-1 # iterates until termination
    global P = flow(P)
    P = [P[x,s,p]/sum(P) for x in 0:size(P,1)-1, s in 0:1, p in 0:1] # renormalization, due to the approximate nature of algorithms.
    P = OffsetArray(P,0:size(P,1)-1,0:1,0:1) # use the offset array
end
# computing the mean and variance of the averaged delay
∑(x) = sum(x)
μ = ∑([(L-x)*P[x,s,p] for x in 1:L, s in 0:1, p in 0:1])
σ² = ∑([(L-x)^2*P[x,s,p] for x in 1:L, s in 0:1, p in 0:1]) - μ^2
μ_c = ∑([s*P[x,s,p] for x in 1:L, s in 0:1, p in 0:1])

if plot_flag == 0
# in order to adapt for multi-threading, for each program, we will use a separate file to store the output.
# if isfile("output.csv")
#     df = CSV.read("output.csv",DataFrame)
#     output = [k_couple,k_uncouple,mean(v_translations),mean(v_transcriptions),mean(v_stalls),mean(v_unstalls),L,ℓ,Eᵦ,E_c, x₀, y₀ ,s₀ , p₀,type, μ, σ², μ_c]
#     push!(df,output)
# else
df = DataFrame(k_couple = [k_couple], k_uncouple = [k_uncouple], v_translations = [v_translation_0], v_transcriptions = [v_transcription_0], v_stalls = [k_stalling_0], v_unstalls = [k_unstalling_0], k_ini_pausings = [k_ini_pausing], L = [L], ℓ = [ℓ], Eᵦ = [Eᵦ], E_c = [E_c], x₀ = [x₀], y₀ = [y₀], s₀ = [s₀], p₀ = [p₀],type=[type], μ = [μ], σ² = [σ²], μ_c = [μ_c])
# end
CSV.write("data/output_$(rand()).csv",df)
using Distributions
t_linspace = [i for i in 0.1:.1:20]
pdf_linspace = [sum([P[x,s,p]*pdf(Erlang(L-x, 1/v_translation_0),t) for x in 1:L for s in 0:1 for p in 0:1]) for t in t_linspace]


# plot(size=(500,500))

label=latexstring("v_{\\textrm{transcription}}=$(v_transcription_0), v_{\\textrm{translation}}=$(v_translation_0), k_{\\textrm{pause}}=$(k_stalling_0), k_{\\textrm{restart}}=$(k_unstalling_0)")
@pgf axis = Axis(
    {
        xlabel = "time/s",
        ylabel = "probability density",
        grid = "major",
        legend_pos  = "north west"
    },
)
x = t_linspace
y = pdf_linspace
# store the data in a dataframe and save it to a csv file
df = DataFrame(t_linspace = x, pdf_linspace = y)
CSV.write("fig/delay/pdf_$(v_transcription_0)_$(v_translation_0)_$(k_stalling_0)_$(k_unstalling_0)_$(k_couple)_$(k_uncouple)_$(E_c)_$(Eᵦ).csv",df)
else
## from distance distribution to waiting time distribution, assuming v_translation is uniform
using Distributions
t_linspace = [i for i in .1:.1:20]
pdf_linspace = [sum([P[x,s,p]*pdf(Erlang(L-x, 1/v_translation_0),t) for x in 1:L for s in 0:1 for p in 0:1]) for t in t_linspace]


# plot(size=(500,500))

label=latexstring("v_{\\textrm{transcription}}=$(v_transcription_0), v_{\\textrm{translation}}=$(v_translation_0), k_{\\textrm{pause}}=$(k_stalling_0), k_{\\textrm{restart}}=$(k_unstalling_0)")
@pgf axis = Axis(
    {
        xlabel = "time/s",
        ylabel = "probability density",
        grid = "major",
        legend_pos  = "north west"
    },
)
x = t_linspace
y = pdf_linspace
# store the data in a dataframe and save it to a csv file
df = DataFrame(t_linspace = x, pdf_linspace = y)
CSV.write("fig/delay/pdf_$(v_transcription_0)_$(v_translation_0)_$(k_stalling_0)_$(k_unstalling_0)_$(k_couple)_$(k_uncouple)_$(E_c)_$(Eᵦ).csv",df)
plot_line = Plot(Coordinates(x,y))
legend_entry = LegendEntry(label)
push!(axis,plot_line)
push!(axis,legend_entry)
# plot!(t_linspace,pdf_linspace,label=label,lw=2)
# xlabel!("time/s")
# ylabel!("probability density")
# title!("waiting time distribution")
# ylims!(-0.002,0.12)

pgfsave("fig/delay/$(v_transcription_0)_$(v_translation_0)_$(k_stalling_0)_$(k_unstalling_0)_$(k_couple)_$(k_uncouple)_$(E_c).svg",axis)
pgfsave("fig/delay/$(v_transcription_0)_$(v_translation_0)_$(k_stalling_0)_$(k_unstalling_0)_$(k_couple)_$(k_uncouple)_$(E_c).tex",axis)
end

