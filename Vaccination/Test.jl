function f(du,u,p,t)
  du[1] = u[2]
  du[2] = -p
end
function condition(u,t,integrator) # Event when event_f(u,t) == 0
  #t ∈ [1:1:10;]
  t>=1 && t==floor(t)
end
function affect!(integrator)
  integrator.u[2] = -integrator.u[2]
end
cb = DiscreteCallback(condition,affect!)

u0 = [50.0,0.0]
tspan = (0.0,15.0)
p = 9.8
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob,Tsit5(),callback=cb,tstops=[1:1:10;])
plot(sol)
