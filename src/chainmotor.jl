"""
Use this structure to simulate chain and motor dynamics. Stiff system, use just Chain to have
faster simulations.
"""
struct ChainMotor{T} <: AbstractChain
    i::Array{T, 1} # Motor states [i, omega]. i is current and omega is rotational speed
    q::Array{T, 1} #generalised coordinates. Contains [x0, alpha_0, alpha_1, alpha_n], the pos of the cart and all link-angles
    u::Array{T, 1} #time derivative of q
    massmat::Array{T, 2}
    fvec::Array{T, 1}
end

function ChainMotor(T::Type)
    i = zeros(T, 1)  # Motor states [i, omega]. i is current and omega is rotational speed
    q = zeros(T, 21) #cart position and all angles are zero
    u = zeros(T, 21) #all velocities are zero
    massmat = zeros(T, 21, 21) #mass matrix, memory is overwritten at each forward simulation step
    fvec = zeros(T, 21) #force vector, memory is overwritten at each forward simulation step
    return ChainMotor{T}(i, q, u, massmat, fvec)
end

ChainMotor() = ChainMotor(Float64)

# Compute state derivatives, with state being [i, omega, x, alpha1, ..., alpha20, dx, dalpha1, ..., dalpha20]
function dxdt_chainmotor(x::AbstractVector{T}, u::T) where {T}
    chain = ChainMotor(T)
    return dxdt_chain!(x, u, chain)
end

function dxdt_chain!(x::AbstractVector, u::T, chain::ChainMotor{T}) where {T}
    chain.i[1] = x[1]
    chain.q .= x[2:22]
    chain.u .= x[23:43]
    return dxdt_chain!(chain, u)
end

function dxdt_chain!(chain::ChainMotor{T}, u::T) where {T}
    # Constant parameters
    k1 = T( 0.263 )  # motor constant
    k2 = T( 0.268 )  # motor constant
    i = T( 3 )  # gear ratio
    r = T( 0.0239 )  #
    R = T( 1.35 )  # armature resistance
    L = T( 1.35*10^-3 )  # armature inductance
    mcart = T( 29.691 )
    mchain = T( 2.3*1.178/20 )
    Ja = T( 24.13*10^-5 )
    Jred = T( Ja + (r/i)^2*( mcart + mchain ))  # reduced moment of inertia

    # Use symmetry of the mass matrix
    symmassmat = Symmetric(chain.massmat, :U)

    # Motor to chain
    Mm = k2*chain.i[1] # k2*i is motor moment
    Fout = i/r*Mm # force by motor on the chain

    # Chain
    quf = vcat(chain.q, chain.u, Fout)
    massmat20!(chain.massmat, quf)
    forcevec20!(chain.fvec, quf)
    dq = chain.u
    du = symmassmat\chain.fvec

    omega = chain.u[1]*i/r
    di::T = -R/L * chain.i[1] - k1/L * omega + u/L

    return vcat(di, dq, du)
end


"""
Defines the dynamics as a differential algebraic equation, as can be solved by DifferentialEquations.jl
"""
function chainmotorDAE!(out::AbstractVector, dx::AbstractVector, x::AbstractVector{T}, u::Number) where {T}
    chain = ChainMotor(T)
    return chainDAE!(out, dx, x, u, chain)
end

function chainDAE!(out::AbstractVector, dx::AbstractVector, x::AbstractVector{T}, u::Number, chain::ChainMotor{T}) where {T}
    chain.i[1] = x[1]
    chain.q .= x[2:22]
    chain.u .= x[23:43]
    return chainDAE!(out, dx, u, chain)
end

function chainDAE!(out::AbstractVector, dx::AbstractVector, u::Number, chain::ChainMotor{T}) where {T}

    # Constant parameters
    k1::T = 0.263# motor constant
    k2::T = 0.268# motor constant
    i::T = 3# gear ratio
    r::T = 0.0239 #
    R::T = 1.35# armature resistance
    L::T = 1.35*10^-3# armature inductance
    mcart::T = 29.691
    mchain::T = 2.3*1.178/20
    Ja::T = 24.13*10^-5
    Jred::T = Ja + (r/i)^2*( mcart + mchain )# reduced moment of inertia

    # Motor to chain
    Mm = k2*chain.i[1] # k2*i is motor moment
    Fout = i/r*Mm # force by motor on the chain

    # Chain
    quf = vcat(chain.q, chain.u, Fout)
    massmat20!(chain.massmat, quf)
    symmassmat = Symmetric(chain.massmat, :U)
    forcevec20!(chain.fvec, quf)
    dq = chain.u

    omega = chain.u[1]*i/r
    di = -R/L * chain.i[1] - k1/L * omega + u/L

    # For DAE formulation
    out[1] = dx[1] - di
    out[2:22] .= dx[2:22] .- dq
    out[23:43] .= symmassmat*dx[23:43] .- chain.fvec
    return out
end
