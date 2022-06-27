mode = "plain"
label = "long_interaction_range" # 
N = 1000 # repeats
k = 30 # v_transcription
ps = [i for i in 2:25] # v_translation
L = 335;
ℓ = 40;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3;
E_c = 1;
k_couple = 100;
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
type = "kinetic_push"
if length(ARGS) == 3 
    E_c = parse(Float64, ARGS[3])
end
label = "long_interaction_range_$(E_c)" #
include("../ttc_plot_ps.jl")