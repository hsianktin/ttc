# type I ttc model, simulation
using Plots
simu_label = "no_barrier"
# simulation of transcription-translation coupling
L = 2400 # length of the sequence
N = 1000 # number of repetitions
ℓ = 100 # maximal distance between the RNAP and the ribosome
v_transcription = 42 # transcription rate nt/second
v_translation = 42 # translation rate nt/second
Eᵦ = 0 # energy cost of uncoupling
x₀ = 0 # initial position of the ribosome
y₀ = 100 # initial position of the RNAP
Γ = []
for i in 1:N
    global x,y = x₀,y₀ # current position of the ribosome and RNAP
    γ = [(x,y)]
    while y < L
        t = rand()
        if x == y
            p = 0 # cannot translate
        elseif x == y - ℓ
            p = v_translation/ (v_translation + v_transcription*exp(-Eᵦ)) # probability of translation
        else
            p = v_translation/(v_translation + v_transcription)

        end
        if t < p
            global x = x + 1 # translate
            push!(γ,(x,y))
        else
            global y = y + 1 # transcribe
            push!(γ,(x,y))
        end
    end
    push!(Γ,γ)
end
plot(legend=false)
for i in 1:N
    plot!(Γ[i],alpha=0.2)
end
xlabel!("ribosome")
ylabel!("RNAP")
savefig("fig/ttc_plot_$(simu_label).png")
histogram([L-Γ[i][end][1] for i in 1:N],label="N=$N")
xlabel!("distance to the end of the sequence")
ylabel!("counts")
title!(simu_label)
savefig("fig/ttc_hist_$(simu_label).png")