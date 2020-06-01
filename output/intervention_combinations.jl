Date(2020,4,7) - Date(2020,3,13) # lockdown 7th April 2020
Date(2020,6,2) - Date(2020,3,13) # School open 2nd June 2020
Date(2020,8,14) - Date(2020,3,13) # Schools closed 14th Augus 2020

Date(2020,8,31) - Date(2020,3,13) # Schools  open 31st august 2020
Date(2020,5,16) - Date(2020,3,13) # Regional lockdown end 16th May
Date(2020,10,23) - Date(2020,3,13) # Schools closed 23rd October 2020

Date(2021,1,4) - Date(2020,3,13)  # schools open 4th Jan 2021
Date(2021,4,9) - Date(2020,3,13)  # Schools closed 9th April 2021
Date(2021,5,3) - Date(2020,3,13)  # Schools open 3rd may 2021
Date(2021,8,6) - Date(2020,3,13)  # Schools closed 6th August 2021
Date(2021,8,30) - Date(2020,3,13)  # Schools open 30th August 2021
Date(2021,10,22) - Date(2020,3,13)  # Schools closed 22nd October 2021


Date(2021,12,31) - Date(2020,3,13)  # we now run intil Dec 2021


#Put in the regional lockdown on day 25 (April 7th)
function regional_lockdown_timing(u,t,integrator)
  !integrator.p.lockdown && t > 25.
end
function affect_regional_lockdown!(integrator)
  integrator.p.T = T_regional_lockdown
  integrator.p.lockdown = true
end
cb_regional_lockdown = DiscreteCallback(regional_lockdown_timing,affect_regional_lockdown!)

#End the regional lockdown on day 64 (May 16th)
function regional_lockdown_ending(u,t,integrator)
  integrator.p.lockdown && t > 64.
end
function affect_regional_lockdown_end!(integrator)
  integrator.p.T = T_normal
  integrator.p.lockdown = false
end
cb_regional_lockdown_end = DiscreteCallback(regional_lockdown_ending,affect_regional_lockdown_end!)

regional_lockdown_starts = CallbackSet(cb_regional_lockdown)
regional_lockdown_starts_and_finishes = CallbackSet(cb_regional_lockdown,cb_regional_lockdown_end)

# First 14 days: Settings for the first 14 days
function first_14_days(u,t,integrator)
    integrator.p.before_week_two && t > 0.
end

# After 14 days to first school opening
function after_first_14_days(u,t,integrator)
     integrator.p.before_week_two && t > 14.
end

#Closure and opening of schools
function open_schools_june(u,t,integrator)
  integrator.p.schools_closed && t > 81.
end

function close_schools_august(u,t,integrator)
  !integrator.p.schools_closed && t > 154.
end

function open_schools_august(u,t,integrator)
  integrator.p.schools_closed && t > 171.
end

function close_schools_october(u,t,integrator)
  !integrator.p.schools_closed && t > 224.
end


function open_schools_jan2021(u,t,integrator)
  integrator.p.schools_closed && t > 297.
end

function close_schools_apr2021(u,t,integrator)
  !integrator.p.schools_closed && t > 392.
end

function open_schools_may2021(u,t,integrator)
  integrator.p.schools_closed && t > 416.
end

function close_schools_aug2021(u,t,integrator)
  !integrator.p.schools_closed && t > 511.
end

function open_schools_aug2021(u,t,integrator)
  integrator.p.schools_closed && t > 535.
end

function close_schools_oct2021(u,t,integrator)
  !integrator.p.schools_closed && t > 588.
end

function affect_first_14_days!(integrator)
   integrator.p.M = 0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school # 0.8 suggested in doc bu +-20% allowed, except schools
end

function affect_open_schools_50pct!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*M_Kenya_school  # the 1.1 in M_Kenya_ho account for 20% increased numbers at home as well as social distancing
  integrator.p.schools_closed = false
end

function affect_open_schools_50pct_candidatesonly!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*candidate_M_school  # the 1.1 in M_Kenya_ho account for 20% increased numbers at home as well as social distancing
  integrator.p.schools_closed = false
end

function affect_open_schools_50pct_primaryonly!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*primary_M_school  # the 1.1 in M_Kenya_ho account for 20% increased numbers at home as well as social distancing
  integrator.p.schools_closed = false
end

function affect_open_schools_50pct_secondaryonly!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*secondary_M_school  # the 1.1 in M_Kenya_ho account for 20% increased numbers at home as well as social distancing
  integrator.p.schools_closed = false
end

function affect_open_schools_50pct_tertiaryonly!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*tertiary_M_school  # the 1.1 in M_Kenya_ho account for 20% increased numbers at home as well as social distancing
  integrator.p.schools_closed = false
end

function affect_open_schools_90pct!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*M_Kenya_school
  integrator.p.schools_closed = false
end

function affect_open_schools_90pct_candidatesonly!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*candidate_M_school
  integrator.p.schools_closed = false
end

function affect_open_schools_90pct_primaryonly!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*primary_M_school  # the 1.1 in M_Kenya_ho account for 20% increased numbers at home as well as social distancing
  integrator.p.schools_closed = false
end

function affect_open_schools_90pct_secondaryonly!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*secondary_M_school  # the 1.1 in M_Kenya_ho account for 20% increased numbers at home as well as social distancing
  integrator.p.schools_closed = false
end

function affect_open_schools_90pct_tertiaryonly!(integrator)
  integrator.p.M = M_Kenya_ho .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*tertiary_M_school  # the 1.1 in M_Kenya_ho account for 20% increased numbers at home as well as social distancing
  integrator.p.schools_closed = false
end

function affect_close_schools!(integrator)  # contact distribution applies to all periods when schooled are closed
  integrator.p.M = 1.2*M_Kenya_ho .+ 0.55*M_Kenya_other .+ 0.55*M_Kenya_work  # contacts at home left at 110% of what they were pre-interventions; room for +-20% allowed by COMORT
  integrator.p.schools_closed = true
  integrator.p.before_week_two = false
end

# open school call backs for the 50% scenarios
cb_open_schools_june_50pct    = DiscreteCallback(open_schools_june,affect_open_schools_50pct!)
cb_open_schools_august_50pct  = DiscreteCallback(open_schools_august,affect_open_schools_50pct!)
cb_open_schools_jan2021_50pct = DiscreteCallback(open_schools_jan2021,affect_open_schools_50pct!)
cb_open_schools_may2021_50pct = DiscreteCallback(open_schools_may2021,affect_open_schools_50pct!)
cb_open_schools_aug2021_50pct = DiscreteCallback(open_schools_aug2021,affect_open_schools_50pct!)

cb_open_schools_june_50pct_candidatesonly    = DiscreteCallback(open_schools_june,affect_open_schools_50pct_candidatesonly!)
cb_open_schools_august_50pct_candidatesonly  = DiscreteCallback(open_schools_august,affect_open_schools_50pct_candidatesonly!)
cb_open_schools_jan2021_50pct_candidatesonly = DiscreteCallback(open_schools_jan2021,affect_open_schools_50pct_candidatesonly!)
cb_open_schools_may2021_50pct_candidatesonly = DiscreteCallback(open_schools_may2021,affect_open_schools_50pct_candidatesonly!)
cb_open_schools_aug2021_50pct_candidatesonly = DiscreteCallback(open_schools_aug2021,affect_open_schools_50pct_candidatesonly!)

cb_open_schools_june_50pct_primaryonly    = DiscreteCallback(open_schools_june,affect_open_schools_50pct_primaryonly!)
cb_open_schools_august_50pct_primaryonly  = DiscreteCallback(open_schools_august,affect_open_schools_50pct_primaryonly!)
cb_open_schools_jan2021_50pct_primaryonly = DiscreteCallback(open_schools_jan2021,affect_open_schools_50pct_primaryonly!)
cb_open_schools_may2021_50pct_primaryonly = DiscreteCallback(open_schools_may2021,affect_open_schools_50pct_primaryonly!)
cb_open_schools_aug2021_50pct_primaryonly = DiscreteCallback(open_schools_aug2021,affect_open_schools_50pct_primaryonly!)

cb_open_schools_june_50pct_secondaryonly    = DiscreteCallback(open_schools_june,affect_open_schools_50pct_secondaryonly!)
cb_open_schools_august_50pct_secondaryonly  = DiscreteCallback(open_schools_august,affect_open_schools_50pct_secondaryonly!)
cb_open_schools_jan2021_50pct_secondaryonly = DiscreteCallback(open_schools_jan2021,affect_open_schools_50pct_secondaryonly!)
cb_open_schools_may2021_50pct_secondaryonly = DiscreteCallback(open_schools_may2021,affect_open_schools_50pct_secondaryonly!)
cb_open_schools_aug2021_50pct_secondaryonly = DiscreteCallback(open_schools_aug2021,affect_open_schools_50pct_secondaryonly!)

cb_open_schools_june_50pct_tertiaryonly    = DiscreteCallback(open_schools_june,affect_open_schools_50pct_tertiaryonly!)
cb_open_schools_august_50pct_tertiaryonly  = DiscreteCallback(open_schools_august,affect_open_schools_50pct_tertiaryonly!)
cb_open_schools_jan2021_50pct_tertiaryonly = DiscreteCallback(open_schools_jan2021,affect_open_schools_50pct_tertiaryonly!)
cb_open_schools_may2021_50pct_tertiaryonly = DiscreteCallback(open_schools_may2021,affect_open_schools_50pct_tertiaryonly!)
cb_open_schools_aug2021_50pct_tertiaryonly = DiscreteCallback(open_schools_aug2021,affect_open_schools_50pct_tertiaryonly!)


# open school call backs for the 90% scenarios
cb_open_schools_june_90pct    = DiscreteCallback(open_schools_june,affect_open_schools_90pct!)
cb_open_schools_august_90pct  = DiscreteCallback(open_schools_august,affect_open_schools_90pct!)
cb_open_schools_jan2021_90pct = DiscreteCallback(open_schools_jan2021,affect_open_schools_90pct!)
cb_open_schools_may2021_90pct = DiscreteCallback(open_schools_may2021,affect_open_schools_90pct!)
cb_open_schools_aug2021_90pct = DiscreteCallback(open_schools_aug2021,affect_open_schools_90pct!)

cb_open_schools_june_90pct_candidatesonly    = DiscreteCallback(open_schools_june,affect_open_schools_90pct_candidatesonly!)
cb_open_schools_august_90pct_candidatesonly  = DiscreteCallback(open_schools_august,affect_open_schools_90pct_candidatesonly!)
cb_open_schools_jan2021_90pct_candidatesonly = DiscreteCallback(open_schools_jan2021,affect_open_schools_90pct_candidatesonly!)
cb_open_schools_may2021_90pct_candidatesonly = DiscreteCallback(open_schools_may2021,affect_open_schools_90pct_candidatesonly!)
cb_open_schools_aug2021_90pct_candidatesonly = DiscreteCallback(open_schools_aug2021,affect_open_schools_90pct_candidatesonly!)

cb_open_schools_june_90pct_primaryonly    = DiscreteCallback(open_schools_june,affect_open_schools_90pct_primaryonly!)
cb_open_schools_august_90pct_primaryonly  = DiscreteCallback(open_schools_august,affect_open_schools_90pct_primaryonly!)
cb_open_schools_jan2021_90pct_primaryonly = DiscreteCallback(open_schools_jan2021,affect_open_schools_90pct_primaryonly!)
cb_open_schools_may2021_90pct_primaryonly = DiscreteCallback(open_schools_may2021,affect_open_schools_90pct_primaryonly!)
cb_open_schools_aug2021_90pct_primaryonly = DiscreteCallback(open_schools_aug2021,affect_open_schools_90pct_primaryonly!)

cb_open_schools_june_90pct_secondaryonly    = DiscreteCallback(open_schools_june,affect_open_schools_90pct_secondaryonly!)
cb_open_schools_august_90pct_secondaryonly  = DiscreteCallback(open_schools_august,affect_open_schools_90pct_secondaryonly!)
cb_open_schools_jan2021_90pct_secondaryonly = DiscreteCallback(open_schools_jan2021,affect_open_schools_90pct_secondaryonly!)
cb_open_schools_may2021_90pct_secondaryonly = DiscreteCallback(open_schools_may2021,affect_open_schools_90pct_secondaryonly!)
cb_open_schools_aug2021_90pct_secondaryonly = DiscreteCallback(open_schools_aug2021,affect_open_schools_90pct_secondaryonly!)

cb_open_schools_june_90pct_tertiaryonly    = DiscreteCallback(open_schools_june,affect_open_schools_90pct_tertiaryonly!)
cb_open_schools_august_90pct_tertiaryonly  = DiscreteCallback(open_schools_august,affect_open_schools_90pct_tertiaryonly!)
cb_open_schools_jan2021_90pct_tertiaryonly = DiscreteCallback(open_schools_jan2021,affect_open_schools_90pct_tertiaryonly!)
cb_open_schools_may2021_90pct_tertiaryonly = DiscreteCallback(open_schools_may2021,affect_open_schools_90pct_tertiaryonly!)
cb_open_schools_aug2021_90pct_tertiaryonly = DiscreteCallback(open_schools_aug2021,affect_open_schools_90pct_tertiaryonly!)

# close school call backs
cb_after_first_14_days = DiscreteCallback(after_first_14_days, affect_close_schools!)
cb_close_schools_august = DiscreteCallback(close_schools_august,affect_close_schools!)
cb_close_schools_october = DiscreteCallback(close_schools_october,affect_close_schools!)
cb_close_schools_apr2021 = DiscreteCallback(close_schools_apr2021,affect_close_schools!)
cb_close_schools_aug2021 = DiscreteCallback(close_schools_aug2021,affect_close_schools!)
cb_close_schools_oct2021 = DiscreteCallback(close_schools_oct2021,affect_close_schools!)


measures_schools_open_june_2020_50pct = CallbackSet(cb_after_first_14_days,
                                            cb_regional_lockdown,  # regional lock down implemented and maintained throughout. Edit this call back to end lock downs
                                            cb_open_schools_june_50pct,
                                            cb_close_schools_august,
                                            cb_open_schools_august_50pct,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021_50pct,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021_50pct,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021_50pct,
                                            cb_close_schools_oct2021)

measures_schools_open_june_2020_50pct_candidatesonly = CallbackSet(cb_after_first_14_days,
                                            cb_regional_lockdown,  # regional lock down implemented and maintained throughout. Edit this call back to end lock downs
                                            cb_open_schools_june_50pct_candidatesonly,
                                            cb_close_schools_august,
                                            cb_open_schools_august_50pct_candidatesonly,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021_50pct_candidatesonly,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021_50pct_candidatesonly,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021_50pct_candidatesonly,
                                            cb_close_schools_oct2021)

measures_schools_open_june_2020_50pct_primaryonly = CallbackSet(cb_after_first_14_days,
                                            cb_regional_lockdown,  # regional lock down implemented and maintained throughout. Edit this call back to end lock downs
                                            cb_open_schools_june_50pct_primaryonly,
                                            cb_close_schools_august,
                                            cb_open_schools_august_50pct_primaryonly,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021_50pct_primaryonly,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021_50pct_primaryonly,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021_50pct_primaryonly,
                                            cb_close_schools_oct2021)

measures_schools_open_june_2020_50pct_secondaryonly = CallbackSet(cb_after_first_14_days,
                                        cb_regional_lockdown,  # regional lock down implemented and maintained throughout. Edit this call back to end lock downs
                                        cb_open_schools_june_50pct_secondaryonly,
                                        cb_close_schools_august,
                                        cb_open_schools_august_50pct_secondaryonly,
                                        cb_close_schools_october,
                                        cb_open_schools_jan2021_50pct_secondaryonly,
                                        cb_close_schools_apr2021,
                                        cb_open_schools_may2021_50pct_secondaryonly,
                                        cb_close_schools_aug2021,
                                        cb_open_schools_aug2021_50pct_secondaryonly,
                                        cb_close_schools_oct2021)

measures_schools_open_june_2020_50pct_tertiaryonly = CallbackSet(cb_after_first_14_days,
                                        cb_regional_lockdown,  # regional lock down implemented and maintained throughout. Edit this call back to end lock downs
                                        cb_open_schools_june_50pct_tertiaryonly,
                                        cb_close_schools_august,
                                        cb_open_schools_august_50pct_tertiaryonly,
                                        cb_close_schools_october,
                                        cb_open_schools_jan2021_50pct_tertiaryonly,
                                        cb_close_schools_apr2021,
                                        cb_open_schools_may2021_50pct_tertiaryonly,
                                        cb_close_schools_aug2021,
                                        cb_open_schools_aug2021_50pct_tertiaryonly,
                                        cb_close_schools_oct2021)

measures_schools_open_june_2020_90pct_candidatesonly = CallbackSet(cb_after_first_14_days,
                                            cb_regional_lockdown,  # regional lock down implemented and maintained throughout. Edit this call back to end lock downs
                                            cb_open_schools_june_90pct_candidatesonly,
                                            cb_close_schools_august,
                                            cb_open_schools_august_90pct_candidatesonly,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021_90pct_candidatesonly,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021_90pct_candidatesonly,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021_90pct_candidatesonly,
                                            cb_close_schools_oct2021)

measures_schools_open_june_2020_90pct_primaryonly = CallbackSet(cb_after_first_14_days,
                                            cb_regional_lockdown,  # regional lock down implemented and maintained throughout. Edit this call back to end lock downs
                                            cb_open_schools_june_90pct_primaryonly,
                                            cb_close_schools_august,
                                            cb_open_schools_august_90pct_primaryonly,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021_90pct_primaryonly,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021_90pct_primaryonly,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021_90pct_primaryonly,
                                            cb_close_schools_oct2021)

measures_schools_open_june_2020_90pct_secondaryonly = CallbackSet(cb_after_first_14_days,
                                        cb_regional_lockdown,  # regional lock down implemented and maintained throughout. Edit this call back to end lock downs
                                        cb_open_schools_june_90pct_secondaryonly,
                                        cb_close_schools_august,
                                        cb_open_schools_august_90pct_secondaryonly,
                                        cb_close_schools_october,
                                        cb_open_schools_jan2021_90pct_secondaryonly,
                                        cb_close_schools_apr2021,
                                        cb_open_schools_may2021_90pct_secondaryonly,
                                        cb_close_schools_aug2021,
                                        cb_open_schools_aug2021_90pct_secondaryonly,
                                        cb_close_schools_oct2021)

measures_schools_open_june_2020_90pct_tertiaryonly = CallbackSet(cb_after_first_14_days,
                                        cb_regional_lockdown,  # regional lock down implemented and maintained throughout. Edit this call back to end lock downs
                                        cb_open_schools_june_90pct_tertiaryonly,
                                        cb_close_schools_august,
                                        cb_open_schools_august_90pct_tertiaryonly,
                                        cb_close_schools_october,
                                        cb_open_schools_jan2021_90pct_tertiaryonly,
                                        cb_close_schools_apr2021,
                                        cb_open_schools_may2021_90pct_tertiaryonly,
                                        cb_close_schools_aug2021,
                                        cb_open_schools_aug2021_90pct_tertiaryonly,
                                        cb_close_schools_oct2021)

measures_schools_open_august_2020_50pct = CallbackSet(cb_after_first_14_days,
                                            cb_regional_lockdown,
                                            cb_open_schools_august_50pct,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021_50pct,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021_50pct,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021_50pct,
                                            cb_close_schools_oct2021)
measures_schools_open_august_2020_50pct_candidatesonly = CallbackSet(cb_after_first_14_days,
                                            cb_regional_lockdown,
                                            cb_open_schools_august_50pct_candidatesonly,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021_50pct_candidatesonly,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021_50pct_candidatesonly,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021_50pct_candidatesonly,
                                            cb_close_schools_oct2021)

measures_schools_open_august_2020_90pct = CallbackSet(cb_after_first_14_days,
                                            cb_regional_lockdown,
                                            cb_open_schools_august_90pct,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021_90pct,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021_90pct,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021_90pct,
                                            cb_close_schools_oct2021)
measures_schools_open_august_2020_90pct_candidatesonly = CallbackSet(cb_after_first_14_days,
                                            cb_regional_lockdown,
                                            cb_open_schools_august_90pct_candidatesonly,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021_90pct_candidatesonly,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021_90pct_candidatesonly,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021_90pct_candidatesonly,
                                            cb_close_schools_oct2021)
