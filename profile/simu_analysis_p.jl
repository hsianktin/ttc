mode = "plain"
N = 1000 # 10s of repeats. The total number is 100.
k = 30 # v_transcription
ps = [i for i in 2:0.2:29] # v_translation
L = 10000;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3;
E_c = 2;
k_couple = 100;
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
d = 27;
type = "kinetic_push"
if length(ARGS) == 3
    E_c = parse(Float64, ARGS[3])
end
label = "analysis_p" # 
# approx_flag = true
append_flag = false
approx_flag = true

include("../ttc_analysis_ps.jl")