using Random
using DataFrames
using CSV
using OffsetArrays
using Statistics
using PGFPlotsX

## initialize parameters
v_transcription_0 = 30.0
v_translation_0 = 20.0
L = 1000 # length of the sequence
ℓ = 4 # characteristic range of couping interactions
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
k_uncouple = k_couple*exp(-Eᵦ) # rate of uncoupling
type = "kinetic_push"
d = 27
# velocity_landscape = DataFrame(translation = Float32[],transcription = Float32[]);
# stalling_landscape = DataFrame(stalling = Float32[],unstalling = Float32[]);
mode = "plain"
N = 10
k_stalling_0 = 0.4 # initial stalling rate
k_unstalling_0 = 0.3 # initial unstalling rate
k_ini_pausing = k_stalling_0
plot_flag = false
# accept ARGS as parameters
if length(ARGS) == 14
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
    d = parse(Float64,ARGS[13])
    type = ARGS[14]
    mode = "plain"
    N = 1000
elseif length(ARGS) == 15
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
    d = parse(Float64,ARGS[13])
    type = ARGS[14]
    mode = ARGS[15]
    N = 1000
elseif length(ARGS) == 16
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
    d = parse(Float64,ARGS[13])
    type = ARGS[14]
    mode = ARGS[15]
    N = parse(Int,ARGS[16])
else
    print("incorrect number of arguments, entering plotting mode")
    plot_flag = true
end
include("ttc_utils.jl")
## Construction the landscapes manually
    # println("mode = plain")
    # println("velocity_landscape.csv not found, initialize a uniform landscape based on base rates")
    velocity_landscape = DataFrame(translation = [v_translation_0 for i in 0:L],transcription = [v_transcription_0 for i in 1:L+1]);
    velocity_landscape.translation[1] = k_translation_initiation # the first step is nothing but the initiation step
    velocity_landscape.transcription[L+1] = k_transcription_termination # the last step is nothing but the termination step
    # # println("stalling_landscape.csv not found")
    stalling_landscape = DataFrame(stalling = [k_stalling_0 for i in 1:L+1],unstalling = [k_unstalling_0 for i in 1:L+1]);
    # "deterministic" stalling at certain position [[@Hatoum2008]] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2819085/
    stalling_landscape.stalling[10] = k_ini_pausing
    # stalling_landscape.unstalling[10] = 0.01

v_translations = velocity_landscape.translation # translation rate nt/second
# println(maximum(v_translations[2:L])- minimum(v_translations[2:L]))
v_translations = OffsetArray(v_translations[1:L+1],0:L) # reset the array, such that v_translations[x] represents the translation elongation or initiation rate at position x
v_transcriptions = velocity_landscape.transcription[1:L+1] # transcription rate nt/second
# for now we assume that y=L is the end of transcription.
v_transcriptions[end-1:end] .= 0.0
v_stalls = stalling_landscape.stalling # stall rate
v_unstalls = stalling_landscape.unstalling # unstall rate

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
        T_0 = Float64[],
        T_1 = Float64[],
        d = Float64[], 
        mean_distance = Float64[],
        mode = String[],
        );

    for i in 1:N
        global x,y,s,p,t = (x₀,y₀,s₀,p₀,0.0)
        push!(v_transcriptions, 0.0)
        push!(v_stalls, 0.0) 
        push!(v_unstalls, 0.0) # after termination, no further transcription
        T0s = []
        T1s = []
        # track the duration of distances larger than a given value d
        local X = Float64[]
        local Y = Float64[]
        local T = Float64[]
        local P = Int[]
        local C = Int[]
        push!(X,x)
        push!(Y,y)
        push!(T,t)
        push!(C,s)
        push!(P,p)
        s₁ = 0
        while x == 0 && y < L
            global x,y,s,p,t = update(x,y,s,p,t,type)
            push!(X,x)
            push!(Y,y)
            push!(T,t)
            push!(C,s)
            push!(P,p)
        end
        if x == 1
            if y - x > d
                s0_flag = true
                t_start = t
            else
                s0_flag = false
                t_start = t
            end
            T_0 = t
            while y < L 
                global x,y,s,p,t = update(x,y,s,p,t,type)
                if s0_flag
                    if y - x <= d && p == 1
                        push!(T0s,t-t_start)
                        s0_flag = false
                        t_start = t
                    end
                else
                    if y - x > d
                        push!(T1s,t-t_start)
                        s0_flag = true
                        t_start = t
                    end
                end
                push!(X,x)
                push!(Y,y)
                push!(T,t)
                push!(C,s)
                push!(P,p)
            end
            local T_transcription = t
            s₁ = s
        else
            s₁ = s
            local T_transcription = t
            while x == 0
                global x,y,s,p,t = update(x,y,s,p,t,type)
                push!(X,x)
                push!(Y,y)
                push!(T,t)
                push!(C,s)
                push!(P,p)
            end
            T_0 = t
        end

        while x < L
            global x,y,s,p,t = update(x,y,s,p,t,type)
            push!(X,x)
            push!(Y,y)
            push!(T,t)
            push!(C,s)
            push!(P,p)
        end
        local T_translation = t #- T_0

        function exposure_duration(X,Y,T,d)
            # calculate the duration at which (Y-X)> d
            # X,Y,T are arrays of the same length
            # d is the threshold
            # return the duration
            # if no such duration exists, return -1
            # if there are multiple such durations, add them up
            # if there are no X>0, return -1
            exposure_duration = 0
            for i in 2:length(X)
                if  Y[i-1] - X[i-1] > d && Y[i-1] < L
                    exposure_duration += T[i] - T[i-1]
                end
            end
            return exposure_duration
        end

        function uncoupled_exposure_duration(X,Y,T,d,C)
            # calculate the duration at which (Y[i]-X[i])> d && C[i] == 0
            exposure_duration = 0
            for i in 2:length(X) 
                if Y[i-1] - X[i-1] > d && C[i] == 0 && Y[i-1] < L
                    exposure_duration += T[i] - T[i-1]
                end
            end
            return exposure_duration
        end

        function mean_distance(X,Y,T)
            # calculate the mean distance between X and Y weighted by (T-T[i-1])
            # X,Y,T are arrays of the same length
            Dδt = 0
            δt = 0
            for i in 2:length(X)
                if Y[i-1] < L
                    Dδt += (Y[i] - X[i]) * (T[i] - T[i-1])
                    δt += T[i] - T[i-1]
                end
            end
            if δt == 0
                return -1
            end
            return Dδt / δt
        end

        local T_exposure = exposure_duration(X,Y,T,d)
        local T_exposure_uncoupled = uncoupled_exposure_duration(X,Y,T,d,C)
        function report_mean(T1s)
            if length(T1s) == 0
                return 0
            else
                return mean(T1s)
            end
        end
        push!(
            df, 
            [
                k_couple,
                k_uncouple, 
                v_translation_0, 
                v_transcription_0, 
                k_stalling_0, 
                k_unstalling_0,
                k_translation_initiation,
                k_ini_pausing, 
                L, 
                ℓ, 
                Eᵦ, 
                E_c, 
                x₀,
                y₀, 
                s₁, # recording the state of coupling at the end of the simulation
                p₀, 
                type, 
                T_transcription, 
                T_translation, 
                T_exposure, 
                T_exposure_uncoupled, 
                report_mean(T0s),
                report_mean(T1s),
                d, 
                mean_distance(X,Y,T),
                mode
                ]
        )
        # df = DataFrame(k_couple = [k_couple], k_uncouple = [k_uncouple], v_translations = [v_translation_0], v_transcriptions = [v_transcription_0], v_stalls = [k_stalling_0], v_unstalls = [k_unstalling_0], k_ini_pausings = [k_ini_pausing], L = [L], ℓ = [ℓ], Eᵦ = [Eᵦ], E_c = [E_c], x₀ = [x₀], y₀ = [y₀], s₀ = [s₀], p₀ = [p₀],type=[type], T_transcription = [T_transcription], T_translation = [T_translation], T_exposure = [T_exposure], T_exposure_uncoupled = [T_exposure_uncoupled], d = [d])
    end
    CSV.write("data/simulation_$(rand()).csv",df)
