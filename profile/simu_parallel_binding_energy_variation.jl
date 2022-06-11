mode = "plain"
label = "binding_energy_variation" # 
N = 1000 # repeats, deprecated, set to 1 for performance.
k = 30 # v_transcription
ps = [i for i in 2:25] # v_translation
L = 335;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦs = [i for i in 1:5];
E_c = 1;
k_couple = 100;
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
type = "kinetic_push"

include("../ttc_plot_ps_Eb.jl")