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
if length(ARGS) == 1
    profile_name = ARGS[1]
    overwrite = false
    include("profile/parallel_$profile_name.jl")
elseif length(ARGS) == 2
    profile_name = ARGS[1]
    overwrite = parse(Bool, ARGS[2])
    include("profile/parallel_$profile_name.jl")
else
    println("Usage: ttc_simu_parallel_run.jl [profile_name]")
end


