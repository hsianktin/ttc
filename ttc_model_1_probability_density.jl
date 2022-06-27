# for a simple model, we can calculate the probability density explicitly.
using Plots
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
# Ps = Array{Array{Float64,3},1}()
# push!(Ps,P₀)
function v_couple(x,y,s,p,k_couple,ℓ,type = "kinetic_push")
    if type == "kinetic_push"
        if abs(x-y) ≤ ℓ && s == 0
            return k_couple
        else
            return 0
        end
    else # type == "soft_couple"
        if s == 0
            return k_couple * exp(-Eᵦ*(y-x)/ℓ)
        else
            return 0
        end
    end
end

function v_uncouple(x,y,s,p,k_uncouple,ℓ,type = "kinetic_push")
    if type == "kinetic_push"
        if abs(x-y) ≤ ℓ && s == 1
            return k_uncouple
        else
            return 0
        end
    else # type == "soft_couple"
        if s == 1
            return k_uncouple * exp(Eᵦ*(y-x)/ℓ)
        else
            return 0
        end
    end
end

function v_translate(x,y,s,p,v_translations,ℓ,type = "kinetic_push")
    if type == "kinetic_push" && x < y
        return v_translations[x]
    elseif x < y # soft_couple
        return v_translations[x]
    else
        return 0
    end
end

function v_transcribe(x,y,s,p,v_transcriptions,ℓ, type = "kinetic_push")
    v = 0
    if p != 0 # if not in stalled state
        if type == "kinetic_push" && y < x + ℓ
            v = v_transcriptions[y]
        elseif type == "soft_couple" # soft_couple
            v = v_transcriptions[y] * exp(-Eᵦ*(y-x)/ℓ)
        else
            v = 0
        end
    end
    return v
end

function v_stall(x,y,s,p,v_stalls,ℓ,type = "kinetic_push")
    if p == 0
        return v_stalls[y]
    else
        return 0
    end
end

function v_unstall(x,y,s,p, v_unstalls, ℓ,type = "kinetic_push")
    if x == y && p == 1
        return v_unstalls[y]*exp(E_c)
    elseif x < y && p == 1
        return v_unstalls[y]
    else
        return 0
    end
end

function transition_kernel(x,y,s,p,ℓ,type = "kinetic_push") # x: position of ribosome, y: position of RNAP, s: state of coupling
    p = zeros(6,1) # probabilities of couple, uncouple, translate, transcribe, stall, unstall
    v = v_couple(x,y,s,p,k_couple,ℓ,type) + v_uncouple(x,y,s,p,k_uncouple,ℓ,type) + v_translate(x,y,s,p,v_translations,ℓ,type) + v_transcribe(x,y,s,p,v_transcriptions,ℓ,type) + v_stall(x,y,s,p,v_stalls,ℓ,type) + v_unstall(x,y,s,p,v_unstalls,ℓ,type)
    p[1] = v_couple(x,y,s,p,k_couple,ℓ,type)/v
    p[2] = v_uncouple(x,y,s,p,k_uncouple,ℓ,type)/v
    p[3] = v_translate(x,y,s,p,v_translations,ℓ,type)/v
    p[4] = v_transcribe(x,y,s,p,v_transcriptions,ℓ,type)/v
    p[5] = v_stall(x,y,s,p,v_stalls,ℓ,type)/v
    p[6] = v_unstall(x,y,s,p,v_unstalls,ℓ,type)/v
    return p
end

function flow(Pₜ) # evolution operator from t to t+1
    # Pₜ is the initial probability distribution
    y,_,_ = size(Pₜ) # reconstruct position of y from the length of probability distribution
    Pₜ₊₁ = zeros(y+1,2,2) # new vector
    Pₜ₊₁ = OffsetArray(Pₜ₊₁,1:y+1,0:1,0:1)
    # one_step_transition
    while sum(Pₜ) > 1e-3
        counter = 0
        for x = 1:y
            for s = 0:1
                for p = 0:1
                    # print("x,y,s,p = $x,$y,$s,$p \r")
                    tk = transition_kernel(x,y,s,p,ℓ,type)
                    if tk[1] > 0 # couple reaction is possible
                        Pₜ[x,s+1,p] += tk[1]*Pₜ[x,s,p]
                    end
                    if tk[2] > 0 # uncouple reaction is possible
                        Pₜ[x,s-1,p] += tk[2]*Pₜ[x,s,p]
                    end
                    if tk[3] > 0 # translate reaction is possible
                        # println("x,s,p=$x,$s,$p")
                        # println(Pₜ[x+1,s,p])
                        # println(Pₜ[x,s,p])
                        Pₜ[x+1,s,p] += tk[3]*Pₜ[x,s,p]
                    end
                    if tk[4] > 0 # transcribe reaction is possible
                        Pₜ₊₁[x,s,p] += tk[4]*Pₜ[x,s,p]
                    end
                    if tk[5] > 0 # stall reaction is possible
                        Pₜ[x,s,p+1] += tk[5]*Pₜ[x,s,p]
                    end
                    if tk[6] > 0 # unstall reaction is possible
                        Pₜ[x,s,p-1] += tk[6]*Pₜ[x,s,p]
                    end
                    # clear the probability of the current position
                    Pₜ[x,s,p] = 0
                end
            end
        end
        counter += 1
    end
    return Pₜ₊₁
end

global P = P₀
@ProgressMeter for i = y₀ + 1 : L
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

