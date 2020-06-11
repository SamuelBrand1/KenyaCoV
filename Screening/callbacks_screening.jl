using LinearAlgebra

############ Utility functions
######### functions for all interventions
function normalize_!(A)
    for neg_i in findall(x->x<0, A)
        A[neg_i]=0
    end
    if sum(A)>0
        return LinearAlgebra.normalize!(A,1)
    end
    A[1]=1
    return A
end
#=function normalize_!(A)
    for neg_i in findall(x->x<0, A)
        A[neg_i]=0
    end
    if size(findall(x->x<0, A),1)==size(A,1)
        A[1]=1
    end
    LinearAlgebra.normalize!(A,1)
    return A
end=#

########### functions for CT
function make_selection_P!(integrator)  # for inplace calculations of probability of an individual of state s to be selected (for testing or as a contact) each day, for a given location and age (∑selection_P[r,a,:] = 1)
    @unpack selection_Pa#=selection per age (for MS)=#,selection_P#=selection per state (for CTH)=#,δ,rel_detection_rate= integrator.p
    u=@view integrator.sol(integrator.t-1)[:,:,:]
    #calculations of Ta and selection_Pa are for MS, while calculations of T and selection_P are for CTH
    for r=1:KenyaCoV_screening.n
        T=sum(u[r,:,1:6])+sum(u[r,:,8])
        for a=1:KenyaCoV_screening.n_a
            Ta=sum(u[r,a,1:6])+u[r,a,8]
            selection_Pa[Int(integrator.t-1),r,a]=Ta/T
            for s=1:8#=n_s=#
                if s==7
                    selection_P[Int(integrator.t-1),r,a,s]=0
                else
                    selection_P[Int(integrator.t-1),r,a,s]=u[r,a,s]/Ta
                end
            end
            normalize_!(selection_P[Int(integrator.t-1),r,a,:])
            if (size(findall(x->x<0, selection_P[Int(integrator.t-1),r,a,:]),1)>0)
                for neg_i in findall(x->x<0, selection_P[Int(integrator.t-1),r,a,:])
                    selection_P[Int(integrator.t-1),r,a,neg_i]=0
                end
                normalize_!(selection_P[Int(integrator.t-1),r,a,:])
            end
        end
        normalize_!(selection_Pa[Int(integrator.t-1),r,:])
    end
end

function make_state_mvt_P!(integrator)#,r,a,prev_day,next_day) make_state_mvt_P
    @unpack state_mvt_P,δ,rel_detection_rate = integrator.p
    prev_day=Int(integrator.t-2);next_day=Int(integrator.t-1)
    sol_next=@view integrator.sol(next_day)[:,:,:];sol_prev=@view integrator.sol(prev_day)[:,:,:]
    for r=1:KenyaCoV_screening.n,a=1:KenyaCoV_screening.n_a
        #calculating numbers of movements between states:
        state_mvt_P[next_day,r,a,1]#=n_r_r=#=sol_prev[r,a,8]      #R(t-1)->R(t)=R(t-1)
        state_mvt_P[next_day,r,a,2]#=KenyaCoV_screening.n_a_r=#=rand(Binomial(Int64(max(sol_next[r,a,8]-sol_prev[r,a,8],0)),1-δ*rel_detection_rate[a]))
        state_mvt_P[next_day,r,a,3]#=n_m_r=#=max(0,sol_next[r,a,8]-sol_prev[r,a,8] - state_mvt_P[next_day,r,a,2]#=KenyaCoV_screening.n_a_r=#)
        state_mvt_P[next_day,r,a,4]#=KenyaCoV_screening.n_a_a=#=max(0,sol_prev[r,a,4] - state_mvt_P[next_day,r,a,2]#=KenyaCoV_screening.n_a_r=#)
        state_mvt_P[next_day,r,a,5]#=n_m_m=#=max(0,sol_prev[r,a,5] - state_mvt_P[next_day,r,a,3]#=n_m_r=#)
        state_mvt_P[next_day,r,a,6]#=n_v_h=#=sol_next[r,a,7] - sol_prev[r,a,7]
        state_mvt_P[next_day,r,a,7]#=n_v_v=#=max(0,sol_prev[r,a,6] - state_mvt_P[next_day,r,a,6]#=n_v_h=#)
        state_mvt_P[next_day,r,a,8]#=n_p_a=#=max(0,sol_next[r,a,4] - sol_prev[r,a,4] + state_mvt_P[next_day,r,a,2]#=KenyaCoV_screening.n_a_r=#)
        state_mvt_P[next_day,r,a,9]#=n_p_m=#=max(0,sol_next[r,a,5] - sol_prev[r,a,5] + state_mvt_P[next_day,r,a,3]#=n_m_r=#)
        state_mvt_P[next_day,r,a,10]#=n_p_v=#=max(0,sol_next[r,a,6] - sol_prev[r,a,6] + state_mvt_P[next_day,r,a,6]#=n_v_h=#)
        state_mvt_P[next_day,r,a,11]#=n_p_p=#=max(0,sol_prev[r,a,3] - state_mvt_P[next_day,r,a,8]#=n_p_a=# - state_mvt_P[next_day,r,a,9]#=n_p_m=# - state_mvt_P[next_day,r,a,10]#=n_p_v=#)
        state_mvt_P[next_day,r,a,12]#=n_e_p=#=max(0,sol_next[r,a,3] - sol_prev[r,a,3] + state_mvt_P[next_day,r,a,8]#=n_p_a=# + state_mvt_P[next_day,r,a,9]#=n_p_m=# + state_mvt_P[next_day,r,a,10]#=n_p_v=#)
        state_mvt_P[next_day,r,a,13]#=n_e_e=#=max(0,sol_prev[r,a,2] - state_mvt_P[next_day,r,a,12]#=n_e_p=#)
        state_mvt_P[next_day,r,a,14]#=n_s_e=#=max(0,sol_next[r,a,2] - sol_prev[r,a,2] + state_mvt_P[next_day,r,a,12]#=n_e_p=#)     #n_s_e2=sol_prev[r,a,1] - sol_next[r,a,1]
        state_mvt_P[next_day,r,a,15]#=n_s_s=#=max(0,sol_prev[r,a,1] - state_mvt_P[next_day,r,a,14]#=n_s_e=#)

        #calculating probability of movements between states:
        if sol_prev[r,a,8]>0   state_mvt_P[next_day,r,a,1]/=sol_prev[r,a,8]    else  state_mvt_P[next_day,r,a,1]=0;end    #p_r_r
        if sol_prev[r,a,4]>0   state_mvt_P[next_day,r,a,2]/=sol_prev[r,a,4]    else  state_mvt_P[next_day,r,a,2]=0;end    #p_a_r
        if sol_prev[r,a,5]>0   state_mvt_P[next_day,r,a,3]/=sol_prev[r,a,5]    else  state_mvt_P[next_day,r,a,3]=0;end    #p_m_r
        if sol_prev[r,a,4]>0   state_mvt_P[next_day,r,a,4]/=sol_prev[r,a,4]    else  state_mvt_P[next_day,r,a,4]=0;end    #p_a_a
        if sol_prev[r,a,5]>0   state_mvt_P[next_day,r,a,5]/=sol_prev[r,a,5]    else  state_mvt_P[next_day,r,a,5]=0;end    #p_m_m
        if sol_prev[r,a,6]>0   state_mvt_P[next_day,r,a,6]/=sol_prev[r,a,6]    else  state_mvt_P[next_day,r,a,6]=0;end    #p_v_h
        if sol_prev[r,a,6]>0   state_mvt_P[next_day,r,a,7]/=sol_prev[r,a,6]    else  state_mvt_P[next_day,r,a,7]=0;end    #p_v_v
        if sol_prev[r,a,3]>0   state_mvt_P[next_day,r,a,8]/=sol_prev[r,a,3]    else  state_mvt_P[next_day,r,a,8]=0;end    #p_p_a
        if sol_prev[r,a,3]>0   state_mvt_P[next_day,r,a,9]/=sol_prev[r,a,3]    else  state_mvt_P[next_day,r,a,9]=0;end    #p_p_m
        if sol_prev[r,a,3]>0   state_mvt_P[next_day,r,a,10]/=sol_prev[r,a,3]    else  state_mvt_P[next_day,r,a,10]=0;end    #p_p_v
        if sol_prev[r,a,3]>0   state_mvt_P[next_day,r,a,11]/=sol_prev[r,a,3]    else  state_mvt_P[next_day,r,a,11]=0;end    #p_p_p
        if sol_prev[r,a,2]>0   state_mvt_P[next_day,r,a,12]/=sol_prev[r,a,2]    else  state_mvt_P[next_day,r,a,12]=0;end    #p_e_p
        if sol_prev[r,a,2]>0   state_mvt_P[next_day,r,a,13]/=sol_prev[r,a,2]    else  state_mvt_P[next_day,r,a,13]=0;end    #p_e_e
        if sol_prev[r,a,1]>0   state_mvt_P[next_day,r,a,14]/=sol_prev[r,a,1]    else  state_mvt_P[next_day,r,a,14]=0;end    #p_s_e
        if sol_prev[r,a,1]>0   state_mvt_P[next_day,r,a,15]/=sol_prev[r,a,1]    else  state_mvt_P[next_day,r,a,15]=0;end    #p_s_s

        #if size(findall(x->x==NaN,state_mvt_P),1)>0   println("ERROR");end
    end
end

function update_states!(integrator,r,a,earliest_day) # updates the states of the numbers of individuals given in V
    @unpack selection_N,state_mvt_P,multinomial_vector2,multinomial_vector4 = integrator.p
    for day=earliest_day-1:-1:1
        if sum(selection_N[:,:])>0
            next_day=Int(integrator.t)-day
            prev_day=Int(integrator.t)-day-1
            if prev_day>0
                #updating Susceptible contacts
                multinomial_vector2 .= rand(Multinomial(selection_N[day+1,1], normalize_!([state_mvt_P[prev_day,r,a,15]#=n_s_s=#,state_mvt_P[prev_day,r,a,14]#=n_s_e=#])))
                selection_N[day,1]+=multinomial_vector2[1]    #S remained S
                selection_N[day,2]+=multinomial_vector2[2]    #S became E
                #updating Exposed contacts
                multinomial_vector2=rand(Multinomial(selection_N[day+1,2], normalize_!([state_mvt_P[prev_day,r,a,13]#=n_e_e=#,state_mvt_P[prev_day,r,a,12]#=n_e_p=#])))
                selection_N[day,2]+=multinomial_vector2[1]   #E remained E
                selection_N[day,3]+=multinomial_vector2[2]   #E became P
                #updating Prodormal contacts
                multinomial_vector4=rand(Multinomial(selection_N[day+1,3], normalize_!([state_mvt_P[prev_day,r,a,11]#=n_p_p=#,state_mvt_P[prev_day,r,a,8]#=n_p_a=#,state_mvt_P[prev_day,r,a,9]#=n_p_m=#,state_mvt_P[prev_day,r,a,10]#=n_p_v=#])))
                selection_N[day,3]+=multinomial_vector4[1]   #P remained P
                selection_N[day,4]+=multinomial_vector4[2]   #P became A
                selection_N[day,5]+=multinomial_vector4[3]   #P became M
                selection_N[day,6]+=multinomial_vector4[4]   #P became V
                #updating Asymptomatic contacts
                multinomial_vector2=rand(Multinomial(selection_N[day+1,4], normalize_!([state_mvt_P[prev_day,r,a,4]#=KenyaCoV_screening.n_a_a=#,state_mvt_P[prev_day,r,a,2]#=KenyaCoV_screening.n_a_r=#])))
                selection_N[day,4]+=multinomial_vector2[1]    #A remained A
                selection_N[day,8]+=multinomial_vector2[2]    #A became R
                #updating Mild symptomatic contacts
                multinomial_vector2=rand(Multinomial(selection_N[day+1,5], normalize_!([state_mvt_P[prev_day,r,a,5]#=n_m_m=#,state_mvt_P[prev_day,r,a,3]#=n_m_r=#])))
                selection_N[day,5]+=multinomial_vector2[1]    #M remained M
                selection_N[day,8]+=multinomial_vector2[2]    #M became R
                #updating seVere symptomatic contacts
                multinomial_vector2=rand(Multinomial(selection_N[day+1,6], normalize_!([state_mvt_P[prev_day,r,a,7]#=n_v_v=#,state_mvt_P[prev_day,r,a,6]#=n_v_h=#])))
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
    @unpack selection_Pa_symptomatics#=selection per age=#,selection_P_symptomatics#=selection per state (M or V for symptomatics)=#,δ,rel_detection_rate= integrator.p
    u=@view integrator.sol(integrator.t-1)[:,:,:]
    for r=1:KenyaCoV_screening.n
        T=sum(u[r,:,5:6])   #number of symptomatics (M and V)
        if T!=0
            for a=1:KenyaCoV_screening.n_a
                selection_Pa_symptomatics[Int(integrator.t-1),r,a]=sum(u[r,a,5:6])/T
                if sum(u[r,a,5:6])>0
                    selection_P_symptomatics[Int(integrator.t-1),r,a,5]=u[r,a,5]/sum(u[r,a,5:6])
                    selection_P_symptomatics[Int(integrator.t-1),r,a,6]=u[r,a,6]/sum(u[r,a,5:6])
                end
                normalize_!(selection_P_symptomatics[Int(integrator.t-1),r,a,:])
            end
        end
        normalize_!(selection_Pa_symptomatics[Int(integrator.t-1),r,:])
    end
end

############
############ Callbacks

############ Common condition method for several intervention callbacks
function timing_daily(u,t,integrator)    t%1==0 && t>0   end

############ Mass Screening intervention
function affect_MS!(integrator)
    @unpack selection_Pa,selection_P = integrator.p
    do_screening!(integrator,selection_Pa,selection_P)
end
function do_screening!(integrator,selection_Pa_vector,selection_P_vector)
    @unpack screening_delay,S_strategy,N_a,toQ,selection_N,detected = integrator.p
    t_test=Int(ceil(integrator.t - screening_delay -1)) #time the test was taken
    if t_test>0
        u_test=@view integrator.sol(t_test)[:,:,:] # state at the testing day
        for r=1:KenyaCoV_screening.n   #region
            if t_test>0 && S_strategy[r,t_test]>0 && sum(selection_Pa_vector[t_test,r,:])>0  #number of tests in this region on the testing day (for which the results will come out today)
                fill!(selection_N,0)   #reset cells to zeros before using this for calculations
                N_a .=  rand(Multinomial(S_strategy[r,t_test], selection_Pa_vector[t_test,r,:]))     #N_a is the number of individuals (ie. tested) per age #N_a[a]
                for a=1:KenyaCoV_screening.n_a
                    if sum(selection_P_vector[t_test,r,a,:])>0
                        #Calculate number of tested per state for a specific age and region
                        selection_N[1+screening_delay,:] .= rand(Multinomial(Int(floor(N_a[a])), selection_P_vector[t_test,r,a,:]))   # the distribution of the tested per state at testing day t_test <=> its index in selection_N is of index 1+screening_delay
                        calculate_positive_tests!(selection_N,integrator.p) #keep only those who tested positive while taking into account the sensitivity/specificity of the test
                        update_states!(integrator,r,a,1+screening_delay)    #updates the states of the contacted people on each day, starting from the testing day (idexed 1+screening_delay in selection_N) and cumulating them to the first (ie. yesterday's) column
                        for s=1:8   #quarantine events SEPAMVR->Q/Qᵥ
                            if selection_N[1,s]>0
                                selection_N[1,s]=min(selection_N[1,s],integrator.sol(integrator.t)[r,a,s])#max_change'ing the events
                                toQ[r,a,s]=selection_N[1,s]
                            end
                        end
                        detected[r,a] = sum(toQ[r,a,:])
                        fill!(selection_N,0)   #reset cells to zeros
                    end
                end
            end
        end
    end
end
callback_MS = DiscreteCallback(timing_daily,affect_MS!)

############ Detection of Hospitalized
"""
This function is a DiscreteCallback integrated to the solver for detecting the Hopitalized
"""
function affect_detect_H!(integrator)
    today=Int(ceil(integrator.t)) #time the detection is performed = today
    integrator.p.detected .= integrator.sol(today)[:,:,12] .- integrator.sol(today-1)[:,:,12] #here we detect all those who moved into the H class
end
callback_detect_H = DiscreteCallback(timing_daily,affect_detect_H!)

############ Contact tracing detecteds
"""
This function is a DiscreteCallback integrated to the solver for the Contact tracing of detecteds
"""
function affect_CT!(integrator)
    @unpack CT_dur,CT_n,state_mvt_P,selection_P,selection_N,M,C₁_set,C₁,C₂,N_a,C₄,toQ,detected = integrator.p
    t=Int(ceil(integrator.t)) #time the contact tracing is performed = today
    if !C₁_set
        C₁=M .* CT_dur
        for a=1:17
            if sum(C₁[:,a])>CT_n
                C₁[:,a]=rand(Multinomial(CT_n, normalize_!(C₁[:,a])))
            end
        end
        C₁_set=true
        @pack! integrator.p = C₁_set,C₁
    end
    sol_t=@view integrator.sol(t)[:,:,:]
    if sum(detected)>0
        for r=1:KenyaCoV_screening.n
            if sum(detected[r,:])>0
                C₂ .= transpose(detected[r,:]) .* C₁
                sum!(N_a,C₂)    #N_a is the number of individuals per age #N_a[a]
                C₄[:,1:CT_dur] .= floor.(N_a ./ CT_dur);    #dispatch the number of contacteds over 14 days
                for c_a=1:KenyaCoV_screening.n_a
                    for d=1:CT_dur
                        if t-d>0 && sum(selection_P[t-d,r,c_a,:])>0
                            selection_N[d,:] .= rand(Multinomial(C₄[c_a,d], selection_P[t-d,r,c_a,:]))        # the distribution of the contacteds per age and per state (total of all days)
                        else
                            selection_N[d,:] .= zeros(Int64,8)
                        end
                    end
                    update_states!(integrator,r,c_a,CT_dur)    #updates the states of the contacted people on each day, cumulating them to the first (ie. yesterday's) column
                    #Quarantine contacts:
                    for s=1:8   #quarantine events SEPAMVR->Q/Qᵥ
                        if selection_N[1,s]>0
                            selection_N[1,s]=min(selection_N[1,s],sol_t[r,c_a,s])#max_change'ing the events
                            toQ[r,c_a,s]=selection_N[1,s]
                        end
                    end
                end
            end
        end
    end
    fill!(detected, 0) #deleting all
end
callback_CT = DiscreteCallback(timing_daily,affect_CT!)

function affect_selection_and_mvt_probabilities!(integrator)#Calculate at each day the probability of an individual to be selected
    if integrator.t>1
        make_selection_P!(integrator)
        make_state_mvt_P!(integrator)
    end
end
callback_selection_and_mvt_probabilities = DiscreteCallback(timing_daily,affect_selection_and_mvt_probabilities!)

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
    if integrator.t>1
        make_selection_P_symptomatics!(integrator)
        make_state_mvt_P!(integrator)
    end
end
callback_selection_and_mvt_probabilities_symptomatics = DiscreteCallback(timing_daily,affect_selection_and_mvt_probabilities_symptomatics!)
