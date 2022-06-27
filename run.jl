###############
# simulation #
###############
# arguments = [
#     "bio_basis", 
#     "bio_control",
#     "bio_test",
#     "frequent_pausing",
#     "long_interaction_range",
#     "slow_coupling",
#     "ls",
#     "Ebs",
#     "Ebs_fast_rib",
#     "Ebs_slow_rib",
#     "alphas",
#     "E_cs",
#     "k_couples"
# ]

# cmds = []
# for argument in arguments
#     cmd = `julia ttc_simu_run.jl $argument false 1`
#     push!(cmds,cmd)
# end

# for cmd in cmds
#     println("running task: $cmd")
#     run(cmd)
# end

###############
## Heat map ##
###############

cmds = []

arguments = [
    # "fast_pausing",
    "slow_pausing"
]

for argument in arguments
    cmd = `julia ttc_heatmap.jl $argument`
    push!(cmds,cmd)
end

for cmd in cmds
    println("running task: $cmd")
    run(cmd)
end

###############
## F_T Plots ##
###############

cmds = []

arguments = [
    "fast_pausing",
    "slow_pausing"
]

for argument in arguments
    cmd = `julia ttc_protection_fraction_p.jl $argument`
    push!(cmds,cmd)
end

for cmd in cmds
    println("running task: $cmd")
    run(cmd)
end


###############
# Delay Dist #
###############

cmds = []
p = 12 # v_translation
q = 30 # v_transcription
L = 335;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3;
E_c = 2;
k_couples = [0,100];
V = [1,10]
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 30;
for k_couple in k_couples
    for v in V
        cmd = `julia ttc_dynamics.jl $q $p $L $ℓ $k_translation_initiation $k_transcription_termination $Eᵦ $E_c $k_couple $(k_stalling_0*v) $(k_unstalling_0*v) $k_ini_pausing  1`
        push!(cmds,cmd)
    end
end

for cmd in cmds
    println("running task: $cmd")
    run(cmd)
end

#################
# Parallel Plot #
#################

cmds = []
arguments = [
    "E_bs_var_ps",
    "E_bs_var_ps_long_range"
]

for argument in arguments
    cmd = `julia ttc_simu_parallel_run.jl $argument true`
    push!(cmds,cmd)
end

map(run,cmds)

#####################
# Delay Bifurcation #
#####################

cmds = []
arguments = [
    "ps_coupled",
    "ps_no_coupling",
    "ps_coupled_long_range",
]

for argument in arguments
    cmd = `julia ttc_delay_bifurcation.jl $argument true`
    push!(cmds,cmd)
end

map(run,cmds)

##################
# Alt Simulation #
##################

cmds = []
arguments = [
    "bio_basis",
    "bio_variant"
]

for argument in arguments
    cmd = `julia ttc_simu_run.jl $argument true 2.0`
    push!(cmds,cmd)
end

map(run,cmds)

#################
# Parallel Alt. #
#################

cmds = []
arguments = [
    "variant",
]

for argument in arguments
    cmd = `julia ttc_simu_parallel_run.jl $argument true`
    push!(cmds,cmd)
end

map(run,cmds)

#################
# Figure 2 #
#################
cmds = []
push!(cmds,`julia ttc_simu_run.jl E_cs false 2.0`)
push!(cmds,`julia ttc_simu_parallel_run.jl ps_var_E_cs false`)
map(run,cmds)

#################
# Figure 3 a - c#
################

cmds = []
arguments = [
    "E_bs_var_ps",
]

for argument in arguments
    cmd = `julia ttc_simu_parallel_run.jl $argument false`
    push!(cmds,cmd)
end

map(run,cmds)


#################
# Figure 3 d - f#
################

#################
# Figure 4 a - c#
################
cmds = []
push!(cmds,`julia ttc_simu_parallel_run.jl alphas_var_ps false`)
push!(cmds,`julia ttc_simu_parallel_run.jl variant false`)
run.(cmds)