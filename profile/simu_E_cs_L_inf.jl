mode = "plain"
label = "E_cs_L_inf"
N = 1000 # repeats
k = 30 # v_transcription
p = 20 # v_translation
L = 10000;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3;
E_cs = [i for i in 0:0.2:5];
k_couple = 100;
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
type = "kinetic_push"
include("../ttc_plot_E_cs.jl")