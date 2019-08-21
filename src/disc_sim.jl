using DifferentialEquations

function chain_odestep(xt::AbstractVector{T}, ut::T, dt::T, chain::AbstractChain = Chain(T)) where {T<:Real}
    # Define ode problem
    tspan = (zero(T), dt)
    differential_vars = [true for i = 1:length(xt)]
    dx0 = zeros(T, length(xt))
    dyn1(out, dx, x, p, t) = chainDAE!(out, dx, x, ut, chain)
    prob = DAEProblem(dyn1, dx0, xt, tspan, differential_vars = differential_vars)

    # solve the problem
    sol = solve(prob, ImplicitMidpoint())

    #return last state
    return sol.u[end]
end
