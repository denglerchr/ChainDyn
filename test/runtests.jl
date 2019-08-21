using ChainDyn

dtype = Float64
chain = Chain(dtype)
chainmotor = ChainMotor(dtype)

xt = zeros(dtype, 42)
xtmotor = zeros(dtype, 43)
ut = dtype(20.0)
dt = dtype(0.01)

xt2 = chain_odestep(xt, ut, dt, chain)
xt2motor = chain_odestep(xtmotor, ut, dt, chainmotor)
