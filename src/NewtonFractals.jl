module NewtonFractals



using Polynomials
using Images

function newtonraphson(x0, f, df, tol=1e-8, maxiter=50)
    
    for i in 1:maxiter
        y = f(x0)
        dy = df(x0)
        
        dx = -y / dy
        x0 = x0 + dx
        if abs(y) < tol
            return x0, i, true
        end
    end
    
    return x0, maxiter, false
    
end


struct NewtonFrac{T <: Real}
    P::Poly{Complex{T}}
    D::Poly{Complex{T}}
    z::Vector{Complex{T}}
    maxiter::Int
    tol::T
    function NewtonFrac(a::AbstractVector{Complex{T}}, maxiter=100, tol=1e-8) where T
        P = Poly{Complex{T}}(a)
        D = polyder(P)
        z = complex(roots(P))
        new{T}(P, D, z, maxiter, T(tol))
    end
end
NewtonFrac(a::AbstractVector{T}, maxiter=100, tol=1e-8) where {T<:Real} = NewtonFrac(complex(a), maxiter, tol)

import Polynomials.roots
roots(N::NewtonFrac) = N.z

nroots(N::NewtonFrac) = length(N.z)

import Polynomials.poly
import Polynomials.polyder

poly(N::NewtonFrac) = N.P
polyder(N::NewtonFrac) = N.D
maxiter(N::NewtonFrac) = N.maxiter
tol(N::NewtonFrac) = N.tol



function newtonraphson(x0, N::NewtonFrac)
    return newtonraphson(x0, poly(N), polyder(N), tol(N), maxiter(N))
end


function cfdcolors(n)
    x = range(0, 4, length=n)
    
    r = zeros(n)
    g = zeros(n)
    b = zeros(n)
    
    for i in 1:n
        xx = x[i]
        if xx <= 1
            r[i] = 1.0
            g[i] = xx
            b[i] = 0.0
        elseif xx <= 2
            r[i] = 1.0 - (xx-1.0)
            g[i] = 1.0
            b[i] = 0.0
        elseif xx <= 3
            r[i] = 0.0
            g[i] = 1.0
            b[i] = xx-2.0
        elseif xx <= 4.0
            r[i] = 0.0
            g[i] = 1.0 - (xx-3.0)
            b[i] = 1.0
        end
    end
    
    return RGB{Float64}.(r, g, b)
end
        
function conv2rgb(z, niter, zlst, clst, imin, imax, amin=0.6)
    
    nz = length(zlst)
    
    errmin = 1e10
    idx = 1
    
    for i in 1:nz
        err = abs(z-zlst[i])
        if err < errmin
            idx = i
            errmin = err
        end
    end
    
    r = clst[idx].r
    g = clst[idx].g
    b = clst[idx].b
    if imin < niter < imax
        a = (imax - niter) / (imax-imin) * (1.0-amin) + amin
    elseif niter <= imin
        a = 1.0
    else
        a = amin
    end
    
    return RGB{Float64}(r*a, g*a, b*a)
end
        


function newtonfrac(N::NewtonFrac{T}, corners; dpi=200,
            itermin=4, itermax=20, amin=0.7, blacknc=true) where T
    xleft, ybottom, xright, ytop = corners[1], corners[2], corners[3], corners[4]

    Δx = xright - xleft
    Δy = ytop - ybottom
    
    nx = round(Int, Δx * dpi)
    ny = round(Int, Δy * dpi)
    
    x = range(T(xleft), T(xright), length=nx)
    y = range(T(ybottom), T(ytop), length=ny)
    
    nz = nroots(N)
    zlst = roots(N)
    clst = cfdcolors(nz)
    
    img = Matrix{RGB{Float64}}(undef, ny, nx)
    
    for i in 1:nx
        xi = x[i]
        for k in 1:ny
            yk = y[k]
            z, niter, conv = newtonraphson(xi + im*yk, N)
            if !conv && blacknc
                img[k,i] = RGB{Float64}(0,0,0)
            else
                img[k,i] = conv2rgb(z, niter, zlst, clst, 
                    itermin, itermax, amin)
            end
        end
    end
    
    return img
end


function newtonfrac(N::NewtonFrac{T}; dpi=200,
        itermin=4, itermax=20, amin=0.7, blacknc=true) where T
    z = roots(N)
    x = real.(z)
    y = imag.(z)
    
    xmin = minimum(x)
    xmax = maximum(x)
    ymin = minimum(y)
    ymax = maximum(y)
    Δx = xmax - xmin
    Δy = ymax - ymin
    
    if Δx < 0.5Δy
        Δx = 0.5Δy
    elseif Δy < 0.5Δx
        Δy = 0.5Δx
    end
    xleft = xmin - 0.5Δx
    xright = xmax + 0.5Δx

    ybottom = ymin - 0.5Δy
    ytop = ymax + 0.5Δy
    
    return newtonfrac(N, (xleft, ybottom, xright, ytop); dpi=dpi, 
        itermin=itermin, itermax=itermax, amin=amin, blacknc=blacknc)

end



end # module
