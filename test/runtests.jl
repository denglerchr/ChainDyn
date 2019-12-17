using ChainDyn

dtype = Float64
chain = Chain(dtype)
chainmotor = ChainMotor(dtype)

xt = zeros(dtype, 42)
xtmotor = zeros(dtype, 43)
ut = dtype(20.0)
dt = dtype(0.01)

@time xt2 = chain_odestep_ODE(xt, ut, dt, chain)
@time xt2motor = chain_odestep_ODE(xtmotor, ut, dt, chainmotor)

@time xt2 = chain_odestep_DAE(xt, ut, dt, chain)
@time xt2motor = chain_odestep_DAE(xtmotor, ut, dt, chainmotor)
