using ChainDyn, DifferentialEquations

function urec(t::T) where {T}
    tau::T = mod(t, 6)
    if 1<tau<3
        return T(20.0)
    elseif 4<tau
        return -T(20.0)
    else
        return T(0.0)
    end
end
############################### Full dynamics##########################
# Create a chain object
chain = ChainMotor(Float64)

# define as a DAE and solve it
x0 = zeros(43)
dx0 = zeros(43)
out = zeros(43)
differential_vars = [true for i = 1:43]
uconst = 10.0
tspan = (0.0, 10.0)

dyn1(out, dx, x, p, t) = chainDAE!(out, dx, x, urec(t), chain)
prob = DAEProblem(dyn1, dx0, x0, tspan, differential_vars = differential_vars)
sol1 = solve(prob, ImplicitMidpoint(), alg_hints = :stiff)

# Visualize
using ModelVisualisations, Plots
pyplot()

X1 = Array{Float64}(undef, 21, length(sol1.u))
t1 = sol1.t
for i = 1:size(X1, 2)
    X1[:, i] .= sol1.u[i][2:22]
end
visualise(:HeavyChain, X1, t1, "dae_full.mp4")

# define as an ODE and solve it
dyn2(x, p, t) = dxdt_chain!(x, urec(t), chain)
prob = ODEProblem(dyn2, x0, tspan)
sol2 = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)

############################### Reduced dynamics##########################
# Create a chain object
chain = Chain(Float64)

# define as a DAE and solve it
x0 = zeros(42)
dx0 = zeros(42)
out = zeros(42)
differential_vars = [true for i = 1:42]
tspan = (0.0, 10.0)

dyn1(out, dx, x, p, t) = chainDAE!(out, dx, x, urec(t), chain)
prob = DAEProblem(dyn1, dx0, x0, tspan, differential_vars = differential_vars)
sol1 = solve(prob, ImplicitMidpoint())

# Visualize
using ModelVisualisations, Plots
pyplot()

X1 = Array{Float64}(undef, 21, length(sol1.u))
t1 = sol1.t
for i = 1:size(X1, 2)
    X1[:, i] .= sol1.u[i][1:21]
end
visualise(:HeavyChain, X1, t1, "dae_red.mp4")

# define as an ODE and solve it
dyn2(x, p, t) = dxdt_chain!(x, urec(t), chain)
prob = ODEProblem(dyn2, x0, tspan)
sol2 = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
