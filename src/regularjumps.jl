"""
For this code to work the 3-dim representation should
be unpacked as a vector
"""

function calculate_infection_rates!(u,p::CoVParameters)
    I_urb_A = @view u[((3-1)*n + 1):((3-1)*n + n)]
    I_urb_D = @view u[((4-1)*n + 1):((4-1)*n + n)]
    I_rur_A = @view u[((3-1)*n + 9*n + 1):((3-1)*n + 9*n + n)]
    I_rur_D = @view u[((4-1)*n + 1):((4-1)*n + n)]
    mul!(p.Î,p.T,I_urb_A .+ I_urb_D  )
    p.Î .+=  I_rur_A .+ I_rur_D
    mul!(p.λ_urb,p.T',p.β .*(p.Î ./N̂))
    p.λ_rur .= p.β .*(p.Î ./N̂)
    return nothing
end

function rates(out,u,p,t)
    calculate_infection_rates!(u,p)
    S_urb = @view u[((1-1)*n + 1):((1-1)*n + n)]
    S_rur = @view u[((1-1)*n + 9*n + 1):((1-1)*n + 9*n + n)]
    E_urb = @view u[((2-1)*n + 1):((2-1)*n + n)]
    E_rur = @view u[((2-1)*n + 9*n + 1):((2-1)*n + 9*n + n)]
    I_urb_A = @view u[((3-1)*n + 1):((3-1)*n + n)]
    I_rur_A = @view u[((3-1)*n + 9*n + 1):((3-1)*n + 9*n + n)]
    I_urb_D = @view u[((4-1)*n + 1):((4-1)*n + n)]
    I_rur_D = @view u[((4-1)*n + 1):((4-1)*n + n)]

    out[1] = (0.1/1000.0)*u[1]*u[2]
    out[2] = 0.01u[2]
end
