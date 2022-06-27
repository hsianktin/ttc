mode = "plain"
N = 1000 # repeats
k = 30 # v_transcription
p = 20 # v_translation
L = 335;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3;
E_c = 1;
k_couples = [10.0^i for i in 0:0.2:2];
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
type = "kinetic_push"
if length(ARGS) == 3 
    E_c = parse(Float64, ARGS[3])
end
label = "k_couples_$(E_c)"
include("../ttc_plot_k_couples.jl")