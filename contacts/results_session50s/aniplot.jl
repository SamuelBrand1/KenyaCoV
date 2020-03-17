using Plots

@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    linewidth --> range(0, 10, length = n)
    alpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end

n = 150
t = range(0, 2π, length = n)
x = sin.(t)
y = cos.(t)

anim = @animate for i ∈ 1:n
    circleplot(x, y, i)
end
gif(anim, "anim_fps15.gif", fps = 15)


function heatgif(A)
    p = heatmap(zeros(size(A, 1),size(A, 2)))
    anim = @animate for i=1:size(A,3)
        heatmap!(p[1], A[:,:,i])
    end
    return anim
end
anim = heatgif(rand(10,10,10))
gif(anim, "./contacts/anim_fps15.gif", fps = 15)
