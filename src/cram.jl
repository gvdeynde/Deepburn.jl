using FFTW
using GenericFFT
using LinearAlgebra
using GenericLinearAlgebra
using Polynomials

include("brent.jl")

"""
Struct that holds a Chebyshev Rational Approximation
"""
struct cra{T <: AbstractFloat}
    name::String
    order::Integer
    rinf::T
    alpha::Vector{Complex{T}}
    theta::Vector{Complex{T}}
    meta::Dict
end

cra(name, order, rinf, alpha, theta) = cra(name, order, rinf, alpha, theta, Dict())

function Base.show(io::IO, c::cra)
    println(io, c.name, "  order = $(c.order)")
    println(io, "r_inf = $(c.rinf)")
    println(io, "alpha")
    for i=1:c.order
        println(io, c.alpha[i])
    end

    println(io, "theta")
    for i=1:c.order
        println(io, c.theta[i])
    end

    println(io, "meta")
    for (key, value) in c.meta
        println(key, ": ", value)
    end
end

"""
Create a polynomial from roots
"""
function fromroots(x)
    n = length(x)
    y = zeros(eltype(x), n+1)
    y[1] = one(eltype(x))
    for j =1:n
        y[2:(j+1)] -= x[j] .* y[1:j]
    end
    return y
end

"""
Create a Hankel matrix based on column vector
"""
function hankel(c)
    n = length(c)
    h = zeros(eltype(c),n,n)
    for i =1:n
        h[1:(n-i+1),i] = c[i:n]
    end
    return h
end

function _CF_Float64(order::Int, start_order::Int,  K::Int, nf::Int)

    size = order-start_order+1

    CFs = Vector{cra}(undef, size)

    one = 1.0
    two = 2.0
    twoj = 2.0im

    w = exp.(twoj*pi.*collect(0:(nf-1))/nf)
    t = real(w)
    scl = 9.0
    F = exp.(scl*((t.-one)./(t.+one)))
    F[isnan.(F)] .= 0.0

    c = real(fft(real(F)))/nf

    f = evalpoly.(w, Ref(c[1:K+1]))

    H = hankel(c[2:K+1])

    SVD = svd(H)

    for n=start_order:order
        s = SVD.S[n+1]
        u = SVD.U[K:-1:1, n+1]
        v = SVD.V[:,n+1]

        zz = zeros(nf-K)
        b = fft(vcat(u,zz))./fft(vcat(v, zz))
        rt = f-s*(w.^K).*b
        zr = roots(Polynomial(v[end:-1:1]))

        qj = zr[findall(x->abs(x)>1,zr)]
        qc = real(fromroots(qj))

        pt = rt.*evalpoly.(w, Ref(qc[end:-1:1]))
        ptc = real(fft(pt)/nf)
        ptc = ptc[n+1:-1:1]

        ci = 0*qj
        for k = 1:n
            q = qj[k]
            q2 = fromroots(qj[findall(x->x!=q,qj)])
            ci[k] = evalpoly(q,ptc[end:-1:1])/evalpoly(q, q2[end:-1:1])
        end

        zi=scl*((qj.-one).^two)./((qj.+one).^two)
        ci=two*two*ci.*zi./(qj.^two .- one)

        rinf = real(one/two*(one+sum(ci./zi)))

        sorder = sortperm(imag(zi))
        zi = zi[sorder]
        ci = ci[sorder]


        CFs[n-start_order+1] = cra("Carathéodory-Fejér", n, rinf, ci, zi, Dict("Type" => "Float64", "K"=>K, "nf"=>nf))

    end

    return CFs

end

function _CF_BigFloat(order::Int, start_order::Int, K::Int, nf::Int, prec::Int)
    oldprec = precision(BigFloat)
    setprecision(BigFloat, prec)

    T = Complex{BigFloat}

    size = order-start_order+1

    CFs = Vector{cra}(undef, size)

    one = T(big"1", big"0")
    two = T(big"2", big"0")
    twoj = T(big"0", big"2")

    bnf = big(nf)*big"1.0"

    w = exp.(twoj*big(pi).*collect(0:big(nf-1))/bnf)
    t = real(w)
    scl = big"9"
    F = exp.(scl*((t.-one)./(t.+one)))
    F[isnan.(F)] .= big"0.0"

    c = real(GenericFFT.generic_fft(real(F)))/bnf

    f = evalpoly.(w, Ref(c[1:K+1]))

    H = hankel(c[2:K+1])

    SVD = svd(H)

    for n=start_order:order
        s = SVD.S[n+1]
        u = SVD.U[K:-1:1, n+1]
        v = SVD.V[:,n+1]

        zz = zeros(nf-K)
        b = GenericFFT.generic_fft(vcat(u,zz))./GenericFFT.generic_fft(vcat(v, zz))
        rt = f-s*(w.^K).*b
        zr = roots(Polynomial(v[end:-1:1]))

        qj = zr[findall(x->abs(x)>1,zr)]
        qc = real(fromroots(qj))

        pt = rt.*evalpoly.(w, Ref(qc[end:-1:1]))
        ptc = real(GenericFFT.generic_fft(pt)/nf)
        ptc = ptc[n+1:-1:1]

        ci = 0*qj
        for k = 1:n
            q = qj[k]
            q2 = fromroots(qj[findall(x->x!=q,qj)])
            ci[k] = evalpoly(q,ptc[end:-1:1])/evalpoly(q, q2[end:-1:1])
        end

        zi=scl*((qj.-one).^big"2")./((qj.+one).^big"2")
        ci=big"4"*ci.*zi./(qj.^big"2" .- one)

        rinf = real(one/two*(one+sum(ci./zi)))

        sorder = sortperm(imag(zi))
        zi = zi[sorder]
        ci = ci[sorder]

        CFs[n-start_order+1] = cra("Carathéodory-Fejér", n, rinf, ci, zi, Dict("Type" => "BigFloat", "K"=>K, "nf"=>nf, "precision"=>precision))

    end

    setprecision(oldprec)

    return CFs
end

function CarathéodoryFejér(order::Int=2; prec::Int=0, start_order::Int=order,  K::Int=75, nf::Int=1024)

    if prec > 0
        res = _CF_BigFloat(order, start_order, K, nf, prec)
    else
        res = _CF_Float64(order, start_order, K, nf)
    end

    if start_order == order
        return res[1]
    else
        return res
    end

end

function (c::cra)(x::Union{Float64, BigFloat, Vector{Float64}, Vector{BigFloat}}, prec::Int=0)

    if prec > 0
        oldprec = precision(BigFloat)
        setprecision(BigFloat, prec)

        two = big"2.0"
        alpha = BigFloat.(string.(c.alpha))
        theta = BigFloat.(string.(c.theta))
        rinf = BigFloat(string(c.rinf))
        xx = BigFloat.(string.(x))

    else
        two = 2.0
        alpha = float.(c.alpha)
        theta = float.(c.theta)
        rinf = float.(c.rinf)
        xx = x
    end

    res = zero(xx)
    n = c.order ÷ 2

    for i=1:n
        res = res .+ alpha[i] ./ (xx .- theta[i])
    end

    res *= two

    if c.order%2 == 1
        res = res .+ alpha[n+1] ./ (xx .- theta[n+1])
    end

    res = res .+ rinf

    if prec > 0
        setprecision(BigFloat, oldprec)
    end

    return real(res)
end

function cra_abs_error(c::cra, x::Union{Float64, BigFloat, Vector{Float64}, Vector{BigFloat}}; prec::Int=0)

    if prec > 0
        oldprec = precision(BigFloat)
        setprecision(BigFloat, prec)

        xx = BigFloat.(string.(x))

    else
        xx = x
    end

    error = exp.(xx) .- c(xx)

    if prec > 0
        setprecision(BigFloat, oldprec)
    end

    return error
end


function cra_rel_error(c::cra, x::Union{Float64, BigFloat, Vector{Float64}, Vector{BigFloat}}; prec::Int=0)

    if prec > 0
        oldprec = precision(BigFloat)
        setprecision(BigFloat, prec)

        xx = BigFloat.(string.(x))

    else
        xx = x
    end

    error = (exp.(xx) .- c(xx)) ./ exp.(xx)

    if prec > 0
        setprecision(BigFloat, oldprec)
    end

    return error
end

function cra_extrema(c::cra; prec::Int=0)
    if prec > 0
        oldprec = precision(BigFloat)
        setprecision(BigFloat, prec)

        extrema = _cra_extrema_Float64(c)

        setprecision(BigFloat, oldprec)
    else
        extrema = _cra_extrema_BigFloat(c)
    end

    return extrema
end

function _cra_extrema_Float64(c::cra)
    
    extrema = zeros(2*c.order)
    
    #Brute force approach
    sz = 4000

    xinit = -(10.) .^ range(pmin,pmax, sz)
    fevals = cra_abs_error(CF, xinit)
    signchange = fevals[1:sz-1].*fevals[2:sz]

    idx = findall(<(0), signchange)

    nrts = length(idx)

    roots = zeros(nrts)

    for i=1:nrts
        roots[i] = brent(t->cra_abs_error(CF,t), xinit[idx[i]], xinit[idx[i]+1])
    end
    sort!(roots)
    
    
    # Find all extrema between two zeros
    for i=1:2*c.order
        opt = optimize(t->-(cra_abs_error(CF,t))^2, roots[i], roots[i+1], Brent(); rel_tol=1e-12)
        extrema[i] = opt.minimizer
    end
    
    return extrema, roots
end
