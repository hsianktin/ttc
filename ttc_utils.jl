## specifying common basic functions

# x:= ribosome position m
# y:= RNAP position n
function v_couple(x,y,s,p,k_couple,ℓ,type = "kinetic_push")
    # if type == "kinetic_push"
        if abs(x-y) ≤ ℓ && s == 0 && x > 0 # x == 0 <=> ribosome not present
            return k_couple
        else
            return 0
        end
    # else # type == "soft_couple"
    #     if s == 0
    #         return k_couple * exp(-Eᵦ*(y-x)/ℓ)
    #     else
    #         return 0
    #     end
    # end
end

function v_uncouple(x,y,s,p,k_uncouple,ℓ,type = "kinetic_push")
    # if type == "kinetic_push"
        if  s == 1
            return k_uncouple
        else
            # if s == 1
            #     println("cannot uncouple, check ($x,$y,$s,$p)")
            # end
            return 0
        end
    # else # type == "soft_couple"
    #     if s == 1
    #         return k_uncouple * exp(Eᵦ*(y-x)/ℓ)
    #     else
    #         return 0
    #     end
    # end
end

function v_translate(x,y,s,p,v_translations,ℓ,type = "kinetic_push")
    # translation should steady compared with transcription
    if type == "kinetic_push" && x < y 
        return v_translations[x]
    elseif x < y # soft_couple. No difference between hard and soft coupling. This way allows for future extension.
        return v_translations[x]
    else
        if p == 0
            # println("collision restricted")
        end
        return 0
    end
end

function v_transcribe(x,y,s,p,v_transcriptions,ℓ, type = "kinetic_push")
    v = 0
    if p == 0 # if not in stalled state
        if s == 1 # in coupled state
            if  y < x + ℓ
                v = v_transcriptions[y]
            else
                if y == x + ℓ
                # println("coupling restricted")
                end
                v = 0
            end
        else # not in coupled state
            v = v_transcriptions[y]
        end
    end
    return v
end

function v_stall(x,y,s,p,v_stalls,ℓ,type = "kinetic_push")
    if p == 0
        return v_stalls[y]
    else
        return 0
    end
end

function v_unstall(x,y,s,p, v_unstalls, ℓ,type = "kinetic_push")
    if type == "kinetic_push"
        if x == y && p == 1 # p(ause) == 1 -> paused, kinetic pushing instead of coupled pushing
            # println("coupled pushing")
            return v_unstalls[y]*exp(E_c)
        elseif x < y && p == 1
            return v_unstalls[y]
        else
            return 0
        end
    elseif type == "chemical_push"
        if x == y && p == 1 & s == 1 # p(ause) == 1 -> paused, chemical pushing instead of coupled pushing
            return v_unstalls[y]*exp(E_c)
        elseif x < y && p == 1
            return v_unstalls[y]
        else
            return 0
        end
    end
end

function transition_kernel(x,y,s,p,ℓ,type = "kinetic_push") # x: position of ribosome, y: position of RNAP, s: state of coupling
    p̄ = zeros(6,1) # probabilities of couple, uncouple, translate, transcribe, stall, unstall
    v = v_couple(x,y,s,p,k_couple,ℓ,type) + v_uncouple(x,y,s,p,k_uncouple,ℓ,type) + v_translate(x,y,s,p,v_translations,ℓ,type) + v_transcribe(x,y,s,p,v_transcriptions,ℓ,type) + v_stall(x,y,s,p,v_stalls,ℓ,type) + v_unstall(x,y,s,p,v_unstalls,ℓ,type)
    p̄[1] = v_couple(x,y,s,p,k_couple,ℓ,type)/v
    p̄[2] = v_uncouple(x,y,s,p,k_uncouple,ℓ,type)/v
    p̄[3] = v_translate(x,y,s,p,v_translations,ℓ,type)/v
    p̄[4] = v_transcribe(x,y,s,p,v_transcriptions,ℓ,type)/v
    p̄[5] = v_stall(x,y,s,p,v_stalls,ℓ,type)/v
    p̄[6] = v_unstall(x,y,s,p,v_unstalls,ℓ,type)/v
    return p̄
end

function update(x,y,s,p,t, type = "kinetic_push")
    p̄ = zeros(6,1) # probabilities of couple, uncouple, translate, transcribe, stall, unstall
    v = v_couple(x,y,s,p,k_couple,ℓ,type) + v_uncouple(x,y,s,p,k_uncouple,ℓ,type) + v_translate(x,y,s,p,v_translations,ℓ,type) + v_transcribe(x,y,s,p,v_transcriptions,ℓ,type) + v_stall(x,y,s,p,v_stalls,ℓ,type) + v_unstall(x,y,s,p,v_unstalls,ℓ,type)
    p̄[1] = v_couple(x,y,s,p,k_couple,ℓ,type)/v
    p̄[2] = v_uncouple(x,y,s,p,k_uncouple,ℓ,type)/v
    p̄[3] = v_translate(x,y,s,p,v_translations,ℓ,type)/v
    p̄[4] = v_transcribe(x,y,s,p,v_transcriptions,ℓ,type)/v
    p̄[5] = v_stall(x,y,s,p,v_stalls,ℓ,type)/v
    p̄[6] = v_unstall(x,y,s,p,v_unstalls,ℓ,type)/v
    # generate waiting time according to the exponential distribution with rate v
    δt = randexp()/v
    # generate a random variable to determine which of the six reaction happens
    r = rand()
    if r < p̄[1]
        # couple
        return x,y,1,p,t+δt
    elseif r < p̄[1]+p̄[2]
        # uncouple
        return x,y,0,p,t+δt
    elseif r < p̄[1]+p̄[2]+p̄[3]
        # translate
        return x+1,y,s,p,t+δt
    elseif r < p̄[1]+p̄[2]+p̄[3]+p̄[4]
        # transcribe
        return x,y+1,s,p,t+δt
    elseif r < p̄[1]+p̄[2]+p̄[3]+p̄[4]+p̄[5]
        # stall
        return x,y,s,1,t+δt
    else
        # unstall
        return x,y,s,0,t+δt
    end
end

function flow(Pₜ) # evolution operator from t to t+1
    # Pₜ is the initial probability distribution
    y,_,_ = size(Pₜ) # reconstruct position of y from the length of probability distribution
    y -= 1 # account for additional "0" state
    Pₜ₊₁ = zeros(y+2,2,2) # new vector
    Pₜ₊₁ = OffsetArray(Pₜ₊₁,0:y+1,0:1,0:1)
    # one_step_transition
    while sum(Pₜ) > 1e-3 # tolerance for inaccuracy.
        counter = 0
        for x = 0:y
            for s = 0:1
                for p = 0:1
                    # print("x,y,s,p = $x,$y,$s,$p \r") # for debugging
                    # print("sum(Pₜ)=$(sum(Pₜ)) \r")
                    tk = transition_kernel(x,y,s,p,ℓ,type)
                    if tk[1] > 0 # couple reaction is possible
                        Pₜ[x,s+1,p] += tk[1]*Pₜ[x,s,p]
                    end
                    if tk[2] > 0 # uncouple reaction is possible
                        Pₜ[x,s-1,p] += tk[2]*Pₜ[x,s,p]
                    end
                    if tk[3] > 0 # translate reaction is possible
                        # println("x,s,p=$x,$s,$p") # for debugging
                        # println(Pₜ[x+1,s,p])
                        # println(Pₜ[x,s,p])
                        Pₜ[x+1,s,p] += tk[3]*Pₜ[x,s,p]
                    end
                    if tk[4] > 0 # transcribe reaction is possible
                        Pₜ₊₁[x,s,p] += tk[4]*Pₜ[x,s,p]
                    end
                    if tk[5] > 0 # stall reaction is possible
                        Pₜ[x,s,p+1] += tk[5]*Pₜ[x,s,p]
                    end
                    if tk[6] > 0 # unstall reaction is possible
                        Pₜ[x,s,p-1] += tk[6]*Pₜ[x,s,p]
                    end
                    # clear the probability of the current position
                    Pₜ[x,s,p] = 0
                end
            end
        end
        counter += 1
    end
    return Pₜ₊₁
end

