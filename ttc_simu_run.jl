using Distributed
using DataFrames
using OffsetArrays
using CSV
using Statistics
using PGFPlotsX
addprocs(7)
begin
    using ProgressMeter
end

## sample parameter profile
label = "bio_control"
N = 10 # repeats
k = 70 # v_transcription
ps = [i for i in 5:70] # v_translation
L = 1000;
ℓ = 2;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 5;
E_c = 30;
k_couple = 100;
k_stalling_0 = 0.3;
k_unstalling_0 = 0.3;
k_ini_pausing = 1e3;
d = 27;
if length(ARGS) == 1
    profile_name = ARGS[1]
    overwrite = false # whether to overwrite existing data
    if isfile("profile/simu_$(profile_name).jl")
        include("profile/simu_$profile_name.jl")
        println("profile $(profile_name) loaded...")
    else
        error("profile $(profile_name) not found...")
    end
elseif length(ARGS) == 2
    profile_name = ARGS[1]
    overwrite = parse(Bool, ARGS[2])
    if isfile("profile/simu_$(profile_name).jl")
        include("profile/simu_$profile_name.jl")
        println("profile $(profile_name) loaded...")
    else
        error("profile $(profile_name) not found...")
    end
elseif length(ARGS) == 3
    profile_name = ARGS[1]
    overwrite = parse(Bool, ARGS[2])
    if isfile("profile/simu_$(profile_name).jl")
        include("profile/simu_$profile_name.jl")
        println("profile $(profile_name) loaded...")
    else
        error("profile $(profile_name) not found...")
    end
else
    println("Usage: ttc_simu_plot.jl [profile_name] [optional: BOOL]")
end