using LinearAlgebra,PoissonRandom

############ Utility functions
######### functions for all interventions
function normalize_!(A)
    for neg_i in findall(x->x<0, A)
        A[neg_i]=0
    end
    if sum(A)>0
        LinearAlgebra.normalize!(A,1)
    else
        A[1]=1
        fill!(A[2:end],0)
    end
    return A
end

########### functions for CT
function make_selection_P!(integrator)  # for inplace calculations of probability of an individual of state s to be selected (for testing or as a contact) each day, for a given location and age (∑selection_P[r,a,:] = 1)
    @unpack selection_Pa#=selection per age (for MS)=#,selection_P#=selection per state (for CTH)=#,dt= integrator.p
    t=Int(integrator.t-dt)
    u=@view integrator.sol(t)[:,:,:]
    #calculations of Ta and selection_Pa are for MS, while calculations of T and selection_P are for CTH
    for r=1:KenyaCoV_screening.n
        T=sum(u[r,:,1:6])+sum(u[r,:,8])
        for a=1:KenyaCoV_screening.n_a
            Ta=sum(u[r,a,1:6])+u[r,a,8]
            selection_Pa[t,r,a]=Ta/T
            for s=1:8#=n_s=#
                if s==7
                    selection_P[t,r,a,s]=0
                else
                    selection_P[t,r,a,s]=u[r,a,s]/Ta
                end
            end
            if (size(findall(x->x<0, selection_P[t,r,a,:]),1)>0)
                for neg_i in findall(x->x<0, selection_P[t,r,a,:])
                    selection_P[t,r,a,neg_i]=0
                end
                #normalize_!(selection_P[t,r,a,:])
            end
            normalize_!(selection_P[t,r,a,:])
        end
        normalize_!(selection_Pa[t,r,:])
    end
end

function make_transition_proba_fb!(integrator)
    @unpack transition_proba_f,transition_proba_b,𝜠,dt = integrator.p
    t⁻=Int(integrator.t-dt-1)
    t⁺=Int(integrator.t-dt)
    if t⁻>0
        u⁻=@view integrator.sol(t⁻)[:,:,1:8]
        u⁺=@view integrator.sol(t⁺)[:,:,1:8]
        for r=1:KenyaCoV_screening.n,a=1:KenyaCoV_screening.n_a
            #calculating forwrds probability of state transition events:
            fill!(transition_proba_f[t⁺,r,a,:],0)
            if u⁻[r,a,1]>0   transition_proba_f[t⁺,r,a,1]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,1]],u⁻[r,a,1])/u⁻[r,a,1],1)    end    #p_s_e
            if u⁻[r,a,2]>0   transition_proba_f[t⁺,r,a,2]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,2]],u⁻[r,a,2])/u⁻[r,a,2],1)    end    #p_e_p
            if u⁻[r,a,3]>0   transition_proba_f[t⁺,r,a,3]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,3]],u⁻[r,a,3])/u⁻[r,a,3],1)    end    #p_p_a
            if u⁻[r,a,3]>0   transition_proba_f[t⁺,r,a,4]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,4]],u⁻[r,a,3])/u⁻[r,a,3],1)    end    #p_p_m
            if u⁻[r,a,3]>0   transition_proba_f[t⁺,r,a,5]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,5]],u⁻[r,a,3])/u⁻[r,a,3],1)    end    #p_p_v
            if u⁻[r,a,6]>0   transition_proba_f[t⁺,r,a,6]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,6]],u⁻[r,a,6])/u⁻[r,a,6],1)    end    #p_v_h
            if u⁻[r,a,5]>0   transition_proba_f[t⁺,r,a,7]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,7]],u⁻[r,a,5])/u⁻[r,a,5],1)    end    #p_m_r
            if u⁻[r,a,4]>0   transition_proba_f[t⁺,r,a,8]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,8]],u⁻[r,a,4])/u⁻[r,a,4],1)    end    #p_a_r

            #calculating backwards probability of state transition events:
            fill!(transition_proba_b[t⁺,r,a,:],0)
            if u⁺[r,a,1]>0   transition_proba_b[t⁺,r,a,1]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,1]],u⁺[r,a,2])/u⁺[r,a,1],1)    end    #p_s_e
            if u⁺[r,a,2]>0   transition_proba_b[t⁺,r,a,2]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,2]],u⁺[r,a,3])/u⁺[r,a,2],1)    end    #p_e_p
            if u⁺[r,a,3]>0   transition_proba_b[t⁺,r,a,3]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,3]],u⁺[r,a,4])/u⁺[r,a,3],1)    end    #p_p_a
            if u⁺[r,a,3]>0   transition_proba_b[t⁺,r,a,4]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,4]],u⁺[r,a,5])/u⁺[r,a,3],1)    end    #p_p_m
            if u⁺[r,a,3]>0   transition_proba_b[t⁺,r,a,5]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,5]],u⁺[r,a,6])/u⁺[r,a,3],1)    end    #p_p_v
            if u⁺[r,a,6]>0   transition_proba_b[t⁺,r,a,6]=min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,6]],u⁺[r,a,7])/u⁺[r,a,6],1)    end    #p_v_h
            if u⁺[r,a,5]>0   transition_proba_b[t⁺,r,a,7]=max(min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,7]],u⁺[r,a,8])/u⁺[r,a,5],1),0)    end    #p_m_r
            if u⁺[r,a,4]>0   transition_proba_b[t⁺,r,a,8]=max(min(min(𝜠[t⁺][KenyaCoV_screening.linear_as_events[r,a,8]],u⁺[r,a,8])/u⁺[r,a,4],1),0)    end    #p_a_r
        end
    end
end

function update_states!(integrator,r,a,earliest_day) # updates the states of the numbers of individuals given in V
    @unpack selection_N,transition_proba_f,multinomial_vector2,multinomial_vector4,dt = integrator.p
    t=Int(integrator.t-dt)
    for day=earliest_day-1:-1:1
        if sum(selection_N[:,:])>0
            next_day=t-day
            prev_day=t-day-1
            if prev_day>0
                #updating Susceptible contacts
                multinomial_vector2 .= rand(Multinomial(selection_N[day+1,1], normalize_!([transition_proba_f[prev_day,r,a,15]#=n_s_s=#,transition_proba_f[prev_day,r,a,14]#=n_s_e=#])))
                selection_N[day,1]+=multinomial_vector2[1]    #S remained S
                selection_N[day,2]+=multinomial_vector2[2]    #S became E
                #updating Exposed contacts
                multinomial_vector2=rand(Multinomial(selection_N[day+1,2], normalize_!([transition_proba_f[prev_day,r,a,13]#=n_e_e=#,transition_proba_f[prev_day,r,a,12]#=n_e_p=#])))
                selection_N[day,2]+=multinomial_vector2[1]   #E remained E
                selection_N[day,3]+=multinomial_vector2[2]   #E became P
                #updating Prodormal contacts
                multinomial_vector4=rand(Multinomial(selection_N[day+1,3], normalize_!([transition_proba_f[prev_day,r,a,11]#=n_p_p=#,transition_proba_f[prev_day,r,a,8]#=n_p_a=#,transition_proba_f[prev_day,r,a,9]#=n_p_m=#,transition_proba_f[prev_day,r,a,10]#=n_p_v=#])))
                selection_N[day,3]+=multinomial_vector4[1]   #P remained P
                selection_N[day,4]+=multinomial_vector4[2]   #P became A
                selection_N[day,5]+=multinomial_vector4[3]   #P became M
                selection_N[day,6]+=multinomial_vector4[4]   #P became V
                #updating Asymptomatic contacts
                multinomial_vector2=rand(Multinomial(selection_N[day+1,4], normalize_!([transition_proba_f[prev_day,r,a,4]#=KenyaCoV_screening.n_a_a=#,transition_proba_f[prev_day,r,a,2]#=KenyaCoV_screening.n_a_r=#])))
                selection_N[day,4]+=multinomial_vector2[1]    #A remained A
                selection_N[day,8]+=multinomial_vector2[2]    #A became R
                #updating Mild symptomatic contacts
                multinomial_vector2=rand(Multinomial(selection_N[day+1,5], normalize_!([transition_proba_f[prev_day,r,a,5]#=n_m_m=#,transition_proba_f[prev_day,r,a,3]#=n_m_r=#])))
                selection_N[day,5]+=multinomial_vector2[1]    #M remained M
                selection_N[day,8]+=multinomial_vector2[2]    #M became R
                #updating seVere symptomatic contacts
                multinomial_vector2=rand(Multinomial(selection_N[day+1,6], normalize_!([transition_proba_f[prev_day,r,a,7]#=n_v_v=#,transition_proba_f[prev_day,r,a,6]#=n_v_h=#])))
                selection_N[day,6]+=multinomial_vector2[1]    #M remained M
                selection_N[day,7]+=multinomial_vector2[2]    #M became R
                #updating hospitalized contacts
                selection_N[day,7]+=selection_N[day+1,7]
                #updating recovered contacts
                selection_N[day,8]+=selection_N[day+1,8]
            end
        end
    end
end

function update_states_contacteds_atLocation!(integrator,r) # updates the states of C₃[r,c_a,s,day]
    @unpack C₃,transition_proba_f,CT_dur,dt = integrator.p
    t=Int(integrator.t-dt)
    for #=r=1:KenyaCoV_screening.n,=#c_a=1:KenyaCoV_screening.n_a,day=2:1:CT_dur
        if sum(C₃[r,c_a,:,day])>0   #If we have some contacteds for this location, age and day
            day⁺=t-CT_dur+day     #day to be updated
            day⁻=t-CT_dur+day-1
            if day⁻>0
                n_s_e=rand(Binomial(C₃[r,c_a,1,day-1],transition_proba_f[day⁻,r,c_a,1]))
                n_e_p=rand(Binomial(C₃[r,c_a,2,day-1],transition_proba_f[day⁻,r,c_a,2]))
                n_p_a=rand(Binomial(C₃[r,c_a,3,day-1],transition_proba_f[day⁻,r,c_a,3]))
                n_p_m=rand(Binomial(C₃[r,c_a,3,day-1],transition_proba_f[day⁻,r,c_a,4]))
                n_p_v=rand(Binomial(C₃[r,c_a,3,day-1],transition_proba_f[day⁻,r,c_a,5]))
                n_v_h=rand(Binomial(C₃[r,c_a,6,day-1],transition_proba_f[day⁻,r,c_a,6]))
                n_m_r=rand(Binomial(C₃[r,c_a,5,day-1],transition_proba_f[day⁻,r,c_a,7]))
                n_a_r=rand(Binomial(C₃[r,c_a,4,day-1],transition_proba_f[day⁻,r,c_a,8]))

                C₃[r,c_a,1,day]+=max(0,C₃[r,c_a,1,day-1]-n_s_e)         #updating Susceptibles
                C₃[r,c_a,2,day]+=max(0,C₃[r,c_a,2,day-1]+n_s_e-n_e_p)   #updating Exposeds
                C₃[r,c_a,3,day]+=max(0,C₃[r,c_a,3,day-1]+n_e_p-n_p_a-n_p_m-n_p_v)   #updating Prodormals
                C₃[r,c_a,4,day]+=max(0,C₃[r,c_a,4,day-1]+n_p_a-n_a_r)   #updating Asymps
                C₃[r,c_a,5,day]+=max(0,C₃[r,c_a,5,day-1]+n_p_m-n_m_r)   #updating Mildly symps
                C₃[r,c_a,6,day]+=max(0,C₃[r,c_a,6,day-1]+n_p_v-n_v_h)   #updating Severe symps
                C₃[r,c_a,7,day]+=max(0,C₃[r,c_a,7,day-1]+n_v_h)         #updating Hospitalized
                C₃[r,c_a,8,day]+=max(0,C₃[r,c_a,8,day-1]+n_a_r+n_m_r)   #updating Recovereds
            end
        end
    end
end
function make_transmission_chain!(integrator)
    @unpack CT_dur,detected_transitions,detected,transition_proba_b,dt = integrator.p
    t=Int(integrator.t-dt)
    #make the history of transmision chain of detecteds using p←
    for r=1:KenyaCoV_screening.n, a=1:KenyaCoV_screening.n_a
        detected_transitions[r,a,CT_dur,:] .= detected[r,a,:]  #get the number of detecteds for all states per location and age
        if sum(detected_transitions[r,a,CT_dur,:])>0   #If we do have some detecteds for this location and age
            for day=CT_dur-1:-1:1   #going backwards in time from day=yesterday=CT_DUR-1 until 14_days_ago=1
                day⁺=t-CT_dur+day+1
                day⁻=t-CT_dur+day     #day to be updated
                if day⁻>0
                    n_s_e=rand(Binomial(detected_transitions[r,a,day+1,2],transition_proba_b[day⁺,r,a,1]))
                    n_e_p=rand(Binomial(detected_transitions[r,a,day+1,3],transition_proba_b[day⁺,r,a,2]))
                    n_p_a=rand(Binomial(detected_transitions[r,a,day+1,4],transition_proba_b[day⁺,r,a,3]))
                    n_p_m=rand(Binomial(detected_transitions[r,a,day+1,5],transition_proba_b[day⁺,r,a,4]))
                    n_p_v=rand(Binomial(detected_transitions[r,a,day+1,6],transition_proba_b[day⁺,r,a,5]))
                    n_v_h=rand(Binomial(detected_transitions[r,a,day+1,7],transition_proba_b[day⁺,r,a,6]))#;print("\ttransition_proba_b[day⁺,r,a,7]=",transition_proba_b[day⁺,r,a,7])
                    n_m_r=rand(Binomial(detected_transitions[r,a,day+1,8],transition_proba_b[day⁺,r,a,7]))
                    n_a_r=rand(Binomial(detected_transitions[r,a,day+1,8],transition_proba_b[day⁺,r,a,8]))

                    detected_transitions[r,a,day,1]+=max(0,detected_transitions[r,a,day+1,1]+n_s_e)         #updating Susceptibles
                    detected_transitions[r,a,day,2]+=max(0,detected_transitions[r,a,day+1,2]-n_s_e+n_e_p)   #updating Exposeds
                    detected_transitions[r,a,day,3]+=max(0,detected_transitions[r,a,day+1,3]-n_e_p+n_p_a+n_p_m+n_p_v)   #updating Prodormals
                    detected_transitions[r,a,day,4]+=max(0,detected_transitions[r,a,day+1,4]-n_p_a+n_a_r)   #updating Asymps
                    detected_transitions[r,a,day,5]+=max(0,detected_transitions[r,a,day+1,5]-n_p_m+n_m_r)   #updating Mildly symps
                    detected_transitions[r,a,day,6]+=max(0,detected_transitions[r,a,day+1,6]-n_p_v+n_v_h)   #updating Severe symps
                    detected_transitions[r,a,day,7]+=max(0,detected_transitions[r,a,day+1,7]-n_v_h)         #updating Hospitalized
                    detected_transitions[r,a,day,8]+=max(0,detected_transitions[r,a,day+1,8]-n_a_r-n_m_r)   #updating Recovereds
                end
            end
        end
    end
    #@pack! integrator.p = detected_transitions
end

######### functions working for MS
function calculate_positive_tests!(V#=matrix [x days,8 states]=#,p) #V represents number of individuals in the different states, this function removes those whose test-result is negative, while taking into account the accuracy and sensitivity of the test
    @unpack test_sensitivity,test_specificity,screening_delay = p
    #keep only those who result in a positive test = {(1-test_specificity)% of {S,R}, test_sensitivity% of EPAMV }
    #all calculations are done for day 1+screening_delay
    V[1+screening_delay,1]=rand(Binomial(V[1+screening_delay,1],1-test_specificity))      #keep false positives from S
    V[1+screening_delay,8]=rand(Binomial(V[1+screening_delay,8],1-test_specificity))      #keep false positives from R
    for s=2:6       #keep true positives from EPAMV
        V[1+screening_delay,s]=rand(Binomial(V[1+screening_delay,s],test_sensitivity))
    end
end

########### functions for Screening Symptomatics
function make_selection_P_symptomatics!(integrator)  # for inplace calculations of probability of an individual of state s to be selected (for testing or as a contact) each day, for a given location and age (∑selection_P[r,a,:] = 1)
    @unpack dt,selection_Pa_symptomatics#=selection per age=#,selection_P_symptomatics#=selection per state (M or V for symptomatics)=#= integrator.p
    t=Int(integrator.t-dt)
    if t>0
        u=@view integrator.sol(t)[:,:,:]
        for r=1:KenyaCoV_screening.n
            T=sum(u[r,:,5:6])   #number of symptomatics (M and V)
            if T!=0
                for a=1:KenyaCoV_screening.n_a
                    selection_Pa_symptomatics[t,r,a]=sum(u[r,a,5:6])/T
                    if sum(u[r,a,5:6])>0
                        selection_P_symptomatics[t,r,a,5]=u[r,a,5]/sum(u[r,a,5:6])
                        selection_P_symptomatics[t,r,a,6]=u[r,a,6]/sum(u[r,a,5:6])
                    end
                    normalize_!(selection_P_symptomatics[t,r,a,:])
                end
            end
            normalize_!(selection_Pa_symptomatics[t,r,:])
        end
    end
end

############
############ Callbacks

############ Common condition method for several intervention callbacks
function timing_daily(u,t,integrator)
    (t-integrator.p.dt)%1==0 && t>=1
end

############ Mass Screening intervention
function affect_MS!(integrator)
    @unpack selection_Pa,selection_P = integrator.p
    do_screening!(integrator,selection_Pa,selection_P)
end
function do_screening!(integrator,selection_Pa_vector,selection_P_vector)
    @unpack screening_delay,strategy,N_a,toQ,selection_N,detected,dt = integrator.p
    t=Int(integrator.t-dt)
    t_test=t - screening_delay #time the test was taken
    if t_test>0
        u_test=@view integrator.sol(t_test)[:,:,:] # state at the testing day
        for r=1:KenyaCoV_screening.n   #region
            if t_test>0 && strategy[r,t_test]>0 && sum(selection_Pa_vector[t_test,r,:])>0  #number of tests in this region on the testing day (for which the results will come out today)
                fill!(selection_N,0)   #reset cells to zeros before using this for calculations
                N_a .=  rand(Multinomial(strategy[r,t_test], LinearAlgebra.normalize!(selection_Pa_vector[t_test,r,:],1)))     #N_a is the number of individuals (ie. tested) per age #N_a[a]
                for a=1:KenyaCoV_screening.n_a
                    if sum(selection_P_vector[t_test,r,a,:])>0
                        #Calculate number of tested per state for a specific age and region
                        selection_N[1+screening_delay,:] .= rand(Multinomial(Int(floor(N_a[a])), LinearAlgebra.normalize!(selection_P_vector[t_test,r,a,:],1)))   # the distribution of the tested per state at testing day t_test <=> its index in selection_N is of index 1+screening_delay
                        calculate_positive_tests!(selection_N,integrator.p) #keep only those who tested positive while taking into account the sensitivity/specificity of the test
                        update_states!(integrator,r,a,1+screening_delay)    #updates the states of the contacted people on each day, starting from the testing day (idexed 1+screening_delay in selection_N) and cumulating them to the first (ie. yesterday's) column
                        for s=1:8   #quarantine events SEPAMVR->Q/Qᵥ
                            if selection_N[1,s]>0
                                selection_N[1,s]=min(selection_N[1,s],integrator.sol(t)[r,a,s])#max_change'ing the events
                                toQ[r,a,s]=max(0,selection_N[1,s])
                                detected[r,a,s]=toQ[r,a,s]
                            end
                        end
                        fill!(selection_N,0)   #reset cells to zeros
                    end
                end
            end
        end
        #fill!(toQ, 0);toQ[28,10,5]=1
        #fill!(detected, 0);detected[28,10,5]=1
        #=if sum(strategy[:,t_test])>0
            print("screened_toQ=",sum(toQ[:,:,:]))
        end=#
    end
end
callback_MS = DiscreteCallback(timing_daily,affect_MS!)

############ Detection of Hospitalized
"""
This function is a DiscreteCallback integrated to the solver for detecting the Hopitalized
"""
function affect_detect_H!(integrator)
    @unpack dt,detected = integrator.p
    today=Int(integrator.t-dt) #time the detection is performed = today
    for r=1:KenyaCoV_screening.n, a=1:KenyaCoV_screening.n_a
        detected[r,a,7] = integrator.sol(today)[r,a,12] - integrator.sol(today-1)[r,a,12] #here we detect all those who moved into the H class
    end
end
callback_detect_H = DiscreteCallback(timing_daily,affect_detect_H!)

############ Contact tracing detecteds
"""
This function is a DiscreteCallback integrated to the solver for the Contact tracing of detecteds
"""
function affect_CT!(integrator)
    @unpack CT_dur,#=CT_n,=#CT_nₜ,selection_P,toQ,detected,detected_transitions,𝜠,𝜠ᵘ,CT_E,β,C₁,C₂,C₃,C₄,Λ,ϵ,ϵ_D,ϵ_V,dt,M_rescaled,M,strategy,test_and_trace_contacts,traced_pop,contacts_to_trace = integrator.p
    t=Int(ceil(integrator.t-dt)) #time the contact tracing is performed = today
    if sum(detected)>0
        #print("\tdeleted before CT=",sum(toQ[:,:,:]))
        #fill!(toQ, 0)                                                                               #####????????????????????????????????????????????????????
        u=@view integrator.sol(t)[:,:,:]
        #make the history of transmision chain of detecteds using p←
        make_transmission_chain!(integrator)
        #Find direct infections to detecteds:
        𝜠ᵘ .= 𝜠[t][KenyaCoV_screening.linear_as_events[:,:,1]]#Set the number of unassigned infections
        #sum_n_traced_contacteds=0;sum_n_actual_contacteds=0;sum_n_actual_contacteds2=0
        for r=1:KenyaCoV_screening.n
            if strategy[r,t]>0
                for a=1:KenyaCoV_screening.n_a,c_a=1:KenyaCoV_screening.n_a,day=1:CT_dur
                    p=0
                    if Λ[t][r,a]>0
                        p=M_rescaled[r,c_a,a]*(ϵ*detected_transitions[r,a,day,3]+ϵ*detected_transitions[r,a,day,4]+ϵ_D*detected_transitions[r,a,day,5]+ϵ_V*detected_transitions[r,a,day,6])/Λ[t][r,a]
                    end
                    nE=rand(Binomial(𝜠ᵘ[r,a],min(1,p)))
                    CT_E[r,c_a,day]+=nE
                    𝜠ᵘ[r,a]-=nE
                end
                #selecting the rest of the contacts (from all states other than the directly exposed)
                                #=
                                CT_E:   actual number of direct exposeds (without pᵗʳᵃᶜᵏ)
                                C₁:     # of actual contacts (i,CT_dur,a,c_a) by all detecteds in i, a
                                C₂:     # of actual contacts (i,c_a,CT_dur) by all detecteds in i and summed for all a
                                C₃:     # of actual contacts (i,c_a,s,CT_dur) distributed by state of contacted
                                        Add CT_E to (i,CT_dur,c_a,2)
                                        Maximize and redistribute
                                        Update contacteds states
                                        Add to isolation schedule
                                =#
                #C₁: # of actual contacts (i,CT_dur,a,c_a) by all detecteds in i, a
                for day=1:CT_dur,a=1:KenyaCoV_screening.n_a
                    C₁[r,day,a,:] .= sum(detected_transitions[r,a,day,:]) .* M[:,a]
                end
                #C₂: # of actual contacts (i,CT_dur,c_a) by all detecteds in i and summed for all a
                for day=1:CT_dur,c_a=1:KenyaCoV_screening.n_a
                    C₂[r,c_a,day] = sum(C₁[r,day,:,c_a])
                end
                #C₃: # of actual contacts (i,CT_dur,c_a,s) distributed by state of contacted
                for c_a=1:KenyaCoV_screening.n_a,day=1:CT_dur
                    if t-CT_dur+day>0
                        C₃[r,c_a,:,day] .= rand(Multinomial(Int(floor(C₂[r,c_a,day])), LinearAlgebra.normalize!(selection_P[t-CT_dur+day,r,c_a,:],1)))
                        C₃[r,c_a,2,day] += CT_E[r,c_a,day]  #Add CT_E to (i,CT_dur,c_a,2)
                    end
                end
                #Rescaling by CT_n per detected
                n_traced_contacteds=pois_rand(CT_nₜ[t]*sum(detected[r,:,:])) #Draw random Poisson of mean CT_n * #detecteds per region (summed for ages, days, and states)
                #print("\t(pois_rand=",n_traced_contacteds,")")
                #sum_n_traced_contacteds+=n_traced_contacteds
                n_actual_contacteds=sum(C₃[r,:,:,:])
                #sum_n_actual_contacteds+=n_actual_contacteds
                if n_actual_contacteds > n_traced_contacteds    #if the number of actual contacteds at location r is > number of traced in r
                    #=for day=1:CT_dur,c_a=1:KenyaCoV_screening.n_a,s=1:8
                        C₃[r,c_a,s,day] = Int(floor(C₃[r,c_a,s,day] * n_traced_contacteds / n_actual_contacteds))
                    end=#
                    C₄ .= reshape(C₃[r,:,:,:],KenyaCoV_screening.n_a*8*CT_dur)
                    C₃[r,:,:,:] .= reshape(  rand(Multinomial(n_traced_contacteds,LinearAlgebra.normalize!(C₄,1)))   ,KenyaCoV_screening.n_a,8,CT_dur)
                end
                #sum_n_actual_contacteds2+=sum(C₃[r,:,:,:])
                update_states_contacteds_atLocation!(integrator,r)    #updates the states of the contacted people on each day, cumulating them to the first (ie. yesterday's) column
                #Quarantine contacts:
                for c_a=1:KenyaCoV_screening.n_a,s=1:8   #quarantine events SEPAMVR->Q/Qᵥ
                    if C₃[r,c_a,s,CT_dur]>0
                        C₃[r,c_a,s,CT_dur]=min(C₃[r,c_a,s,CT_dur],u[r,c_a,s])   #max_change'ing the events
                        toQ[r,c_a,s]+=C₃[r,c_a,s,CT_dur]
                        if test_and_trace_contacts && contacts_to_trace#For when CT of contacts is performed
                            traced_pop[r,c_a,s,1]=C₃[r,c_a,s,CT_dur]
                        end
                    end
                end
            end
        end
        fill!(detected, 0) #deleting all
        fill!(detected_transitions, 0) #deleting all
        #print("\tsum_traced=",sum_n_traced_contacteds,"\tsum_actual=",sum_n_actual_contacteds,"\tsum_actual2=",sum_n_actual_contacteds2)
        #print("\tCT_E=",sum(CT_E))
        fill!(CT_E,0) #deleting all

        #print("\tCT toQ=",sum(toQ[:,:,:]))
    end
    @pack! integrator.p = contacts_to_trace
end
callback_CT = DiscreteCallback(timing_daily,affect_CT!)

function affect_selection_and_mvt_probabilities!(integrator)#Calculate at each day the probability of an individual to be selected
    make_selection_P!(integrator)
    make_transition_proba_fb!(integrator)
end
callback_selection_and_mvt_probabilities = DiscreteCallback(timing_daily,affect_selection_and_mvt_probabilities!)


###########
function affect_CT2!(integrator)
    @unpack CT_dur,#=CT_n,=#CT_nₜ,selection_P,toQ,detected,detected_transitions,𝜠,𝜠ᵘ,CT_E,β,C₁,C₂,C₃,C₄,Λ,ϵ,ϵ_D,ϵ_V,dt,M_rescaled,M,strategy,test_and_trace_contacts,traced_pop,contacts_to_trace = integrator.p
    t=Int(ceil(integrator.t-dt)) #time the contact tracing is performed = today
    if sum(detected)>0
        #print("\tdeleted before CT=",sum(toQ[:,:,:]))
        #fill!(toQ, 0)                                                                               #####????????????????????????????????????????????????????
        u=@view integrator.sol(t)[:,:,:]
        #make the history of transmision chain of detecteds using p←
        make_transmission_chain!(integrator)
        #Find direct infections to detecteds:
        𝜠ᵘ .= 𝜠[t][KenyaCoV_screening.linear_as_events[:,:,1]]#Set the number of unassigned infections
        #sum_n_traced_contacteds=0;sum_n_actual_contacteds=0;sum_n_actual_contacteds2=0
        for r=1:KenyaCoV_screening.n
            if strategy[r,t]>0
                for a=1:KenyaCoV_screening.n_a,c_a=1:KenyaCoV_screening.n_a,day=1:CT_dur
                    p=0
                    if Λ[t][r,a]>0
                        p=M_rescaled[r,c_a,a]*(ϵ*detected_transitions[r,a,day,3]+ϵ*detected_transitions[r,a,day,4]+ϵ_D*detected_transitions[r,a,day,5]+ϵ_V*detected_transitions[r,a,day,6])/Λ[t][r,a]
                    end
                    nE=rand(Binomial(𝜠ᵘ[r,a],min(1,p)))
                    CT_E[r,c_a,day]+=nE
                    𝜠ᵘ[r,a]-=nE
                end
                #selecting the rest of the contacts (from all states other than the directly exposed)
                                #=
                                CT_E:   actual number of direct exposeds (without pᵗʳᵃᶜᵏ)
                                C₁:     # of actual contacts (i,CT_dur,a,c_a) by all detecteds in i, a
                                C₂:     # of actual contacts (i,c_a,CT_dur) by all detecteds in i and summed for all a
                                C₃:     # of actual contacts (i,c_a,s,CT_dur) distributed by state of contacted
                                        Add CT_E to (i,CT_dur,c_a,2)
                                        Maximize and redistribute
                                        Update contacteds states
                                        Add to isolation schedule
                                =#
                #C₁: # of actual contacts (i,CT_dur,a,c_a) by all detecteds in i, a
                for day=1:CT_dur,a=1:KenyaCoV_screening.n_a
                    C₁[r,day,a,:] .= sum(detected_transitions[r,a,day,:]) .* M[:,a]
                end
                #C₂: # of actual contacts (i,CT_dur,c_a) by all detecteds in i and summed for all a
                for day=1:CT_dur,c_a=1:KenyaCoV_screening.n_a
                    C₂[r,c_a,day] = sum(C₁[r,day,:,c_a])
                end
                #C₃: # of actual contacts (i,CT_dur,c_a,s) distributed by state of contacted
                for c_a=1:KenyaCoV_screening.n_a,day=1:CT_dur
                    if t-CT_dur+day>0
                        C₃[r,c_a,:,day] .= rand(Multinomial(Int(floor(C₂[r,c_a,day])), LinearAlgebra.normalize!(selection_P[t-CT_dur+day,r,c_a,:],1)))
                        C₃[r,c_a,2,day] += CT_E[r,c_a,day]  #Add CT_E to (i,CT_dur,c_a,2)
                    end
                end
                #Rescaling by CT_n per detected
                n_traced_contacteds=pois_rand(CT_nₜ[t]*sum(detected[r,:,:])) #Draw random Poisson of mean CT_n * #detecteds per region (summed for ages, days, and states)
                #print("\t(pois_rand=",n_traced_contacteds,")")
                #sum_n_traced_contacteds+=n_traced_contacteds
                n_actual_contacteds=sum(C₃[r,:,:,:])
                #sum_n_actual_contacteds+=n_actual_contacteds
                if n_actual_contacteds > n_traced_contacteds    #if the number of actual contacteds at location r is > number of traced in r
                    #=for day=1:CT_dur,c_a=1:KenyaCoV_screening.n_a,s=1:8
                        C₃[r,c_a,s,day] = Int(floor(C₃[r,c_a,s,day] * n_traced_contacteds / n_actual_contacteds))
                    end=#
                    C₄ .= reshape(C₃[r,:,:,:],KenyaCoV_screening.n_a*8*CT_dur)
                    C₃[r,:,:,:] .= reshape(  rand(Multinomial(n_traced_contacteds,LinearAlgebra.normalize!(C₄,1)))   ,KenyaCoV_screening.n_a,8,CT_dur)
                end
                #sum_n_actual_contacteds2+=sum(C₃[r,:,:,:])
                update_states_contacteds_atLocation!(integrator,r)    #updates the states of the contacted people on each day, cumulating them to the first (ie. yesterday's) column
                #Quarantine contacts:
                for c_a=1:KenyaCoV_screening.n_a,s=1:8   #quarantine events SEPAMVR->Q/Qᵥ
                    if C₃[r,c_a,s,CT_dur]>0
                        C₃[r,c_a,s,CT_dur]=min(C₃[r,c_a,s,CT_dur],u[r,c_a,s])   #max_change'ing the events
                        toQ[r,c_a,s]+=C₃[r,c_a,s,CT_dur]
                        #=if test_and_trace_contacts && contacts_to_trace#For when CT of contacts is performed
                            traced_pop[r,c_a,s,1]=C₃[r,c_a,s,CT_dur]
                        end=#
                    end
                end
            end
        end
        fill!(detected, 0) #deleting all
        fill!(detected_transitions, 0) #deleting all
        #print("\tsum_traced=",sum_n_traced_contacteds,"\tsum_actual=",sum_n_actual_contacteds,"\tsum_actual2=",sum_n_actual_contacteds2)
        #print("\tCT_E=",sum(CT_E))
        fill!(CT_E,0) #deleting all

        #print("\tCT toQ=",sum(toQ[:,:,:]))
    end
    @pack! integrator.p = contacts_to_trace
end
callback_CT2 = DiscreteCallback(timing_daily,affect_CT2!)

############ Detection of Symptomatics
"""
This function is a DiscreteCallback integrated to the solver for screening symptomatics
"""
function affect_SympS!(integrator)
    @unpack selection_Pa_symptomatics,selection_P_symptomatics = integrator.p
    do_screening!(integrator,selection_Pa_symptomatics,selection_P_symptomatics)
end
callback_SympS = DiscreteCallback(timing_daily,affect_SympS!)

function affect_selection_and_mvt_probabilities_symptomatics!(integrator)#Calculate at each day the probability of an individual to be selected
    if integrator.t>0
        make_selection_P_symptomatics!(integrator)
        make_transition_proba_fb!(integrator)
    end
end
callback_selection_and_mvt_probabilities_symptomatics = DiscreteCallback(timing_daily,affect_selection_and_mvt_probabilities_symptomatics!)


########### Contact tracing of traced contacts (CTC)
function affect_ScreenPopulation!(integrator)
    @unpack traced_pop = integrator.p
            ######
            ## traced_pop[n,n_a,8,1+screening_delay]=[r,a,s,day since tested]: contains the contacts that were identified and will be tested
            ##
    ##First, we update the state of those waiting to be tested (in traced_pop[r,a,s,1:end-1])
    update_states_in_traced_pop!(integrator)

    ##Second, those whose test results come out today (in traced_pop[r,a,s,end]) get tested and are moved to detected=traced_pop⁺[r,a,s] (we only store positive testers)
    calculate_positive_tests_from_traced_pop!(integrator)

    ##Third, we contact trace those in detected=traced_pop⁺ by calling callback_CT
end
callback_screen_population = DiscreteCallback(timing_daily,affect_ScreenPopulation!)

function update_states_in_traced_pop!(integrator) # updates the states of traced_pop[r,c_a,s,day]
    @unpack traced_pop,transition_proba_f#=,CT_dur=#,dt,screening_delay = integrator.p
    t=Int(integrator.t-dt) #today
    #print(t,"--BEF-traced_pop[:,:,2:6,:]=",sum(traced_pop[:,:,2:6,:]),"\t--BEF-traced_pop[:,:,2:6,2]=",sum(traced_pop[:,:,2:6,2]))
    for r=1:KenyaCoV_screening.n,c_a=1:KenyaCoV_screening.n_a,day=screening_delay+1:-1:2# screening_delay<=>end
        n_s_e=rand(Binomial(traced_pop[r,c_a,1,day-1],transition_proba_f[t,r,c_a,1]))
        n_e_p=rand(Binomial(traced_pop[r,c_a,2,day-1],transition_proba_f[t,r,c_a,2]))
        n_p_a=rand(Binomial(traced_pop[r,c_a,3,day-1],transition_proba_f[t,r,c_a,3]))
        n_p_m=rand(Binomial(traced_pop[r,c_a,3,day-1],transition_proba_f[t,r,c_a,4]))
        n_p_v=rand(Binomial(traced_pop[r,c_a,3,day-1],transition_proba_f[t,r,c_a,5]))
        n_v_h=rand(Binomial(traced_pop[r,c_a,6,day-1],transition_proba_f[t,r,c_a,6]))
        n_m_r=rand(Binomial(traced_pop[r,c_a,5,day-1],transition_proba_f[t,r,c_a,7]))
        n_a_r=rand(Binomial(traced_pop[r,c_a,4,day-1],transition_proba_f[t,r,c_a,8]))

        traced_pop[r,c_a,1,day]=max(0,traced_pop[r,c_a,1,day-1]-n_s_e)         #updating Susceptibles
        traced_pop[r,c_a,2,day]=max(0,traced_pop[r,c_a,2,day-1]+n_s_e-n_e_p)   #updating Exposeds
        traced_pop[r,c_a,3,day]=max(0,traced_pop[r,c_a,3,day-1]+n_e_p-n_p_a-n_p_m-n_p_v)   #updating Prodormals
        traced_pop[r,c_a,4,day]=max(0,traced_pop[r,c_a,4,day-1]+n_p_a-n_a_r)   #updating Asymps
        traced_pop[r,c_a,5,day]=max(0,traced_pop[r,c_a,5,day-1]+n_p_m-n_m_r)   #updating Mildly symps
        traced_pop[r,c_a,6,day]=max(0,traced_pop[r,c_a,6,day-1]+n_p_v-n_v_h)   #updating Severe symps
        traced_pop[r,c_a,7,day]=max(0,traced_pop[r,c_a,7,day-1]+n_v_h)         #updating Hospitalized
        traced_pop[r,c_a,8,day]=max(0,traced_pop[r,c_a,8,day-1]+n_a_r+n_m_r)   #updating Recovereds
        for s=1:8 # those updated from yesterday to today are deleted from yesterday
            traced_pop[r,c_a,s,day-1]=0
        end
    end
    #print("\tAFT-traced_pop[:,:,2:6,:]=",sum(traced_pop[:,:,2:6,:]),"\t--BEF-traced_pop[:,:,2:6,2]=",sum(traced_pop[:,:,2:6,2]));println()
end

function calculate_positive_tests_from_traced_pop!(integrator)
    @unpack test_sensitivity,test_specificity,screening_delay,traced_pop,detected = integrator.p
    #keep only those who result in a positive test = {(1-test_specificity)% of {S,R}, test_sensitivity% of EPAMV }
    #all calculations are done for day 1+screening_delay
    today=screening_delay+1 #index of today = those whose test results come out today
    #print("traced_pop[:,:,2:6,today]=",sum(traced_pop[:,:,2:6,today]),"\ttraced_pop[:,:,5:6,today]=",sum(traced_pop[:,:,5:6,today]))
    for r=1:KenyaCoV_screening.n,c_a=1:KenyaCoV_screening.n_a#,day=screening_delay+1:-1:2# screening_delay<=>end
        detected[r,c_a,1]=rand(Binomial(traced_pop[r,c_a,1,today],1-test_specificity))      #keep false positives from S
        detected[r,c_a,8]=rand(Binomial(traced_pop[r,c_a,8,today],1-test_specificity))      #keep false positives from S
        for s=2:6   #keep true positives from EPAMV
            detected[r,c_a,s]=rand(Binomial(traced_pop[r,c_a,s,today],test_sensitivity))
        end
        for s=1:8   #those whose test result came out were moved to the detected vector and removed from traced_pop
            traced_pop[r,c_a,s,today]=0
        end
    end
    #print("\tdetected[:,:,2:6]=",sum(detected[:,:,2:6]),"\tdetected[:,:,5:6]=",sum(detected[:,:,5:6]));println()
end
