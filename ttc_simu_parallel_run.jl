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

A = [1,2,3]
push!(A, 4)
using Plots
x = [i for i in 1:5]
y = [i for i in 1:5]
plot(x,y)
plot!(2x,y)

x .+ 1
x .- 1

