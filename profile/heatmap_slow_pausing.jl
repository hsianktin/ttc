mode = "plain"
min = 3
step = 0.5
max = 30
qs = [k for k in min:step:max]
ps = [p for p in min:step:max]
L = 335;
ℓ = 4;
k_translation_initiation = 1.0;
k_transcription_termination = 0.1;
Eᵦ = 3;
E_c = 2;
k_couple = 100;
k_stalling_0 = 0.4;
k_unstalling_0 = 0.3;
k_ini_pausing = k_stalling_0;
label = "slow_pausing_L_$(L)"
N = 1000
d = 27
type = "kinetic_push"