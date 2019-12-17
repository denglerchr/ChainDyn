using ChainDyn, Plots, ModelVisualisations

dtype = Float64
chain = Chain(dtype; friction_link = 1e-4)
dt = dtype(0.01)
allt = 0:dt:10

# Define controller
r(x) = -20*x[1]+5*x[2]

# Start state
x0 = zeros(42)
x0[1] = 0.5

# Simulate
X = zeros(dtype, length(x0), length(allt))
X[:, 1] .= x0
@time for i = 1:length(allt)-1
    obs = vcat( X[1, i] , X[22, i] , X[2, i] , X[23, i] )
    X[:, i+1] .= chain_odestep_ODE(X[:, i], r(obs), dt, chain)
end
pyplot()
visualise(:HeavyChain, Float64.(X[1:21, :]), Float64.(allt), "chain.mp4")
