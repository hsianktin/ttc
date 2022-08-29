mode = "plain"
N = 1000 # repeats
k = 30 # v_transcription
p = 20 # v_translation
L = 1000;
ℓ = 4;
αs = [10.0^(i) for i in -2:0.25:2];
k_transcription_termination = 0.1;
Eᵦ = 3.0;
E_c = 2.0;
k_couple = 100;
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
if length(ARGS) == 3 
    E_c = parse(Float64, ARGS[3])
end
label = "alphas_approx"
type = "kinetic_push"
approx_flag = true
include("../ttc_plot_alphas.jl")