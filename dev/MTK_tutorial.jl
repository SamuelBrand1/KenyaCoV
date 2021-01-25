push!(LOAD_PATH,joinpath(homedir(),"Github/KenyaCoV/src"))
# using Revise
# import KenyaCoV
using DifferentialEquations,ModelingToolkit,Latexify,SparseArrays,StaticArrays,LinearAlgebra,Plots
#Tutorial code for ModelingToolkit

@variables x y
z = x^2 + y
A = [x^2+y 0 2x
     0     0 2y
     y^2+x 0 0]

sparse(A)

function f(u)
    [u[1]-u[3],u[1]^2-u[2],u[3]+u[2]]
end
f([x,y,z])

@variables u[1:3]
@variables(S[1:47,1:17],
            I[1:47,1:17])
f(u)

#Building functions: build_function takes in two args: 1) symbolic expr or array of symbolic expressions
#2) arguments for the expression

to_compute = [x^2 + y, y^2 + x]
f_expr = build_function(to_compute,[x,y])
myf = eval(f_expr[1]) #<----- Returns the computed array
@time myf(SA[2.,3.])

myf! = eval(f_expr[2])
out = zeros(2)
myf!(out,[3.,4.])#<-----Inplace operation
out
write("test.jl", string(f_expr[2]))
f_expr = build_function(to_compute,[x,y],expression=Val{false})

myg = eval(f_expr[1])
@time myg(SA[1.,2.])

#Building Non-Allocating Parallel Functions for Sparse Matrices
@variables x y
N = 10000
A = sparse(Tridiagonal([x^i for i in 1:N-1],[x^i * y^(8-i) for i in 1:N], [y^i for i in 1:N-1]))
T = sparse(Tridiagonal([1. for i in 1:N-1],[1. for i in 1:N], [1. for i in 1:N-1]))
# f_par_expr = build_function(A,[x,y],parallel=ModelingToolkit.MultithreadedForm())
f_expr = build_function(A,[x,y],parallel=ModelingToolkit.MultithreadedForm())[2];
f_expr_np = build_function(A,[x,y])[2];
myf = eval(f_expr)
myf_np = eval(f_expr_np)
@time myf(T,[1.,2.])
@time myf_np(T,[1.,2.])
#Derivatives
@variables t
@derivatives D'~t
z = t + t^2
D(z)
expand_derivatives(D(z))
ModelingToolkit.jacobian([x+x*y,x^2+y],[x,y])

@variables S[1:2] I[1:2]
A = [1. 2.
    2. 1.]
λ = A*I
rate = -S.*λ
f_expr = build_function(rate,S,I)
myf = eval(f_expr[1])
myf([100.,100.],[1.,1.])

ModelingToolkit.jacobian(rate,[S;I])

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(D(x)) ~ (y-x),
       D(y) ~ x*(1-z)-y,
       D(z) ~ x*y - β*z]

sys = ODESystem(eqs)
sys = ode_order_lowering(sys)

u0 = [D(x) => 2.0,
      x => 1.0,
      y => 0.0,
      z => 0.0]

p  = [σ => 28.0,
      ρ => 10.0,
      β => 8/3]

tspan = (0.0,100.0)
prob = ODEProblem(sys,u0,tspan,p,jac=true)
@time sol = solve(prob,Rosenbrock23())
using Plots; plot(sol,vars = (y,z))
