"""
Use this structure to simulate chain dynamics without motor dynamics. Use ChainMotor if you want
to include motor dynamics, but this will turn simulations lower
"""
struct Chain{T}
    q::Array{T, 1} #generalised coordinates. Contains [x0, alpha_0, alpha_1, alpha_n], the pos of the cart and all link-angles
    u::Array{T, 1} #time derivative of q
    massmat::Array{T, 2}
    fvec::Array{T, 1}
end


function Chain(T::Type)
    q = zeros(T, 21) #cart position and all angles are zero
    u = zeros(T, 21) #all velocities are zero
    massmat = zeros(T, 21, 21) #mass matrix, memory is overwritten at each forward simulation step
    fvec = zeros(T, 21) #force vector, memory is overwritten at each forward simulation step
    return Chain{T}(q, u, massmat, fvec)
end

Chain() = Chain(Float64)

# Compute state derivatives, with state being [i, omega, x, alpha1, ..., alpha20, dx, dalpha1, ..., dalpha20]
function dxdt_chain(x::AbstractVector{T}, u::T) where {T}
    chain = Chain(T)
    return dxdt_chain!(x, u, chain)
end

function dxdt_chain!(x::AbstractVector, u::T, chain::Chain{T}) where {T}
    chain.q .= x[1:21]
    chain.u .= x[22:42]
    return dxdt_chain!(chain, u)
end

function dxdt_chain!(chain::Chain{T}, u::T) where {T}
    # Constant parameters
    k1 = T( 0.263 )  # motor constant
    k2 = T( 0.268 )  # motor constant
    i = T( 3 )  # gear ratio
    r = T( 0.0239 )  #
    R = T( 1.35 )  # armature resistance
    mcart = T( 29.691 )
    mchain = T( 2.3*1.178/20 )
    Ja = T( 24.13*10^-5 )
    Jred = T( Ja + (r/i)^2*( mcart + mchain ))  # reduced moment of inertia

    # Motor to chain
    omega = chain.u[1]*i/r # roatational speed of the motor anker
    Mm = k2/R*(u-k1*omega)
    Fout = i/r*Mm # force by motor on the chain

    # Chain
    quf = vcat(chain.q, chain.u, Fout)
    massmat20!(chain.massmat, quf)
    symmassmat = Symmetric(chain.massmat, :U)
    forcevec20!(chain.fvec, quf)
    dq = chain.u
    du = symmassmat\chain.fvec

    return vcat(dq, du)
end


"""
Defines the dynamics as a differential algebraic equation, as can be solved by DifferentialEquations.jl
"""
function chainDAE!(out::AbstractVector, dx::AbstractVector, x::AbstractVector{T}, u::Number) where {T}
    chain = Chain(T)
    return chainDAE!(out, dx, x, u, chain)
end

function chainDAE!(out::AbstractVector, dx::AbstractVector, x::AbstractVector{T}, u::Number, chain::Chain{T}) where {T}
    chain.q .= x[1:21]
    chain.u .= x[22:42]
    return chainDAE!(out, dx, u, chain)
end

function chainDAE!(out::AbstractVector, dx::AbstractVector, u::Number, chain::Chain{T}) where {T}

    # Constant parameters
    k1::T = 0.263# motor constant
    k2::T = 0.268# motor constant
    i::T = 3# gear ratio
    r::T = 0.0239 #
    R::T = 1.35# armature resistance
    mcart::T = 29.691
    mchain::T = 2.3*1.178/20
    Ja::T = 24.13*10^-5
    Jred::T = Ja + (r/i)^2*( mcart + mchain )# reduced moment of inertia

    # Motor to chain
    omega = chain.u[1]*i/r # roatational speed of the motor anker
    Mm = k2/R*(u-k1*omega)
    Fout = i/r*Mm # force by motor on the chain

    # Chain
    quf = vcat(chain.q, chain.u, Fout)
    massmat20!(chain.massmat, quf)
    symmassmat = Symmetric(chain.massmat, :U)
    forcevec20!(chain.fvec, quf)
    dq = chain.u

    # For DAE formulation
    out[1:21] .= dx[1:21] .- dq
    out[22:42] .= symmassmat*dx[22:42] .- chain.fvec
    return out
end
