# This script plots the results from test_flight_interventions.jl

ũ_no_controls = [reshape(u,n,n_s,2) for u in sol_tl_no_controls.u]
cum_I_no_controls =  [sum(u[:,7:8,:]) for u in ũ_no_controls ]
ũ_complete_shutdown = [reshape(u,n,n_s,2) for u in sol_tl_complete_shutdown.u]
cum_I_complete_shutdown =  [sum(u[:,7:8,:]) for u in ũ_complete_shutdown ]
ũ_halved_after_week = [reshape(u,n,n_s,2) for u in sol_tl_halved_after_week.u]
cum_I_halved_after_week =  [sum(u[:,7:8,:]) for u in ũ_halved_after_week ]
ũ_quartered_after_week = [reshape(u,n,n_s,2) for u in sol_tl_quartered_after_week.u]
cum_I_quartered_after_week =  [sum(u[:,7:8,:]) for u in ũ_quartered_after_week ]
ũ_shutdown_after_week = [reshape(u,n,n_s,2) for u in sol_tl_shutdown_after_week.u]
cum_I_shutdown_after_week =  [sum(u[:,7:8,:]) for u in ũ_shutdown_after_week ]

plot(sol_tl_no_controls.t[cum_I_no_controls.>0],cum_I_no_controls[cum_I_no_controls.>0],yaxis=:log,lab="No control measures")
plot!(sol_tl_halved_after_week.t[cum_I_halved_after_week.>0],cum_I_halved_after_week[cum_I_halved_after_week.>0],yaxis=:log,lab="Flights halved after one week")
plot!(sol_tl_quartered_after_week.t[cum_I_quartered_after_week.>0],cum_I_quartered_after_week[cum_I_quartered_after_week.>0],yaxis=:log,lab="Flights quartered after one week")
plot!(sol_tl_shutdown_after_week.t[cum_I_shutdown_after_week.>0],cum_I_shutdown_after_week[cum_I_shutdown_after_week.>0],yaxis=:log,lab="All flights stopped after one week")
plot!(sol_tl_complete_shutdown.t[cum_I_complete_shutdown.>0],cum_I_complete_shutdown[cum_I_complete_shutdown.>0],yaxis=:log,lab="All flights stopped immediately")
xlabel!("Day")
ylabel!("Cumulative incidence")
title!("External prevalence 1e-7")
