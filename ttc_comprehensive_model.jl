# for a simple model, we can calculate the probability density explicitly.
# the changes not merged yet.
using ProgressMeter
using DataFrames
using CSV
using OffsetArrays
using Statistics

## initialize parameters
L = 1000 # length of the sequence
ℓ = 100 # characteristic range of couping interactions
k_translation_initiation = 1.0 # translation initiation rate
k_transcription_termination = 0.1 # transcription termination rate
Eᵦ = 1 # energy cost of uncoupling
E_c = 1 # energy contributing to pushing RNAP out of stalling state
# v_y′ = v_y * exp(-Eᵦ)
x₀ = 1 # initial position of the ribosome
y₀ = 100 # initial position of the RNAP
s₀ = 0 # initially uncoupled state
p₀ = 1 # initially stalled RNAP
k_couple = 100 # rate of coupling
k_uncouple = 100*exp(-Eᵦ) # rate of uncoupling
type = "kinetic_push"
# velocity_landscape = DataFrame(translation = Float32[],transcription = Float32[]);
# stalling_landscape = DataFrame(stalling = Float32[],unstalling = Float32[]);
v_translation_0 = 42 # initial translation velocity
v_transcription_0 = 42 # initial transcription velocity
k_stalling_0 = 0.1 # initial stalling rate
k_unstalling_0 = 0.1 # initial unstalling rate
if isfile("velocity_landscape.csv")
    velocity_landscape = CSV.read("velocity_landscape.csv",DataFrame)
else
    println("velocity_landscape.csv not found, initialize a uniform landscape based on base rates")
    velocity_landscape = DataFrame(translation = [v_translation_0 for i in 1:L],transcription = [v_transcription_0 for i in 1:L]);
    CSV.write("velocity_landscape.csv",velocity_landscape)
end
if isfile("stalling_landscape.csv")
    stalling_landscape = CSV.read("stalling_landscape.csv",DataFrame)
else
    println("stalling_landscape.csv not found")
    stalling_landscape = DataFrame(stalling = [k_stalling_0 for i in 1:L],unstalling = [k_unstalling_0 for i in 1:L]);
    stalling_landscape.stalling[y₀] = 1e2
    stalling_landscape.unstalling[y₀] = 0.1
    CSV.write("stalling_landscape.csv",stalling_landscape)
end
v_translations = velocity_landscape.translation # translation rate nt/second
v_transcriptions = velocity_landscape.transcription # transcription rate nt/second
v_stalls = stalling_landscape.stalling # stall rate
v_unstalls = stalling_landscape.unstalling # unstall rate

## initialize probability distribution assuming known initial x₀ and y₀
P₀ = zeros(Float64,y₀,2,2)
P₀ = OffsetArray(P₀,1:y₀,0:1,0:1)
P₀[x₀,s₀,p₀] = 1
P = P₀

# common basic functions are just included
include("ttc_utils.jl")

global P = P₀
@showprogress for i = y₀ + 1 : L
    global P = flow(P)
    P = [P[x,s,p]/sum(P) for x in 1:size(P,1), s in 0:1, p in 0:1] # renormalization
    P = OffsetArray(P,1:size(P,1),0:1,0:1)
    # push!(Ps,P)
end
# computing the mean and variance of the averaged delay
∑(x) = sum(x)
μ = ∑([(L-x)*P[x,s,p] for x in 1:L, s in 0:1, p in 0:1])
σ² = ∑([(L-x)^2*P[x,s,p] for x in 1:L, s in 0:1, p in 0:1]) - μ^2
μ_c = ∑([s*P[x,s,p] for x in 1:L, s in 0:1, p in 0:1])

if isfile("output.csv")
    df = CSV.read("output.csv",DataFrame)
    output = [k_couple,k_uncouple,mean(v_translations),mean(v_transcriptions),mean(v_stalls),mean(v_unstalls),L,ℓ,Eᵦ,E_c, x₀, y₀ ,s₀ , p₀,type, μ, σ², μ_c]
    push!(df,output)
else
    df = DataFrame(k_couple = [k_couple], k_uncouple = [k_uncouple], v_translations = [mean(v_translations)], v_transcriptions = [mean(v_transcriptions)], v_stalls = [mean(v_stalls)], v_unstalls = [mean(v_unstalls)], L = [L], ℓ = [ℓ], Eᵦ = [Eᵦ], E_c = [E_c], x₀ = [x₀], y₀ = [y₀], s₀ = [s₀], p₀ = [p₀],type=[type], μ = [μ], σ² = [σ²], μ_c = [μ_c])
end

CSV.write("output.csv",df)

