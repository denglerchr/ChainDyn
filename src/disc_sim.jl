using DifferentialEquations, Sundials

function chain_odestep_DAE(xt::AbstractVector{T}, ut::T, dt::T, chain::ChainDyn.AbstractChain = Chain(T)) where {T<:Real}
    # Define ode problem
    tspan = (zero(T), dt)
    differential_vars = [true for i = 1:length(xt)]
    dx0 = zeros(T, length(xt))
    dyn1(out, dx, x, p, t) = chainDAE!(out, dx, x, ut, chain)
    prob = DAEProblem(dyn1, dx0, xt, tspan, differential_vars = differential_vars)

    # solve the problem
    sol = solve(prob, IDA(); alg_hints = :stiff)

    #return last state
    return sol.u[end]
end

function chain_odestep_ODE(xt::AbstractVector{T}, ut::T, dt::T, chain::ChainDyn.AbstractChain = Chain(T); alg = DP5()) where {T<:Real}
    # Define ode problem
    tspan = (zero(T), dt)
    dyn2(x, p, t) = dxdt_chain!(x, ut, chain)
    prob = ODEProblem(dyn2, xt, tspan)
    sol = solve(prob, alg ,reltol=1e-8,abstol=1e-8)

    #return last state
    return sol.u[end]
end
