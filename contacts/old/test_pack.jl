using Parameters
@with_kw mutable struct S2
    a::Int=0
    b=[]
    c=[1,2,3]
end

function f(p)
    @unpack a,b,c = p
    for i=1:10  a+=1  end
    push!(b,a)
    c[2]=a
    @pack! p=a
end

p=S2()
f(p)
println(p.c)
