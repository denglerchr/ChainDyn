module ChainDyn

using LinearAlgebra

abstract type AbstractChain end # general chain type

include("friction.jl")

include("forcevec20.jl")
include("massmat20.jl")

include("chainmotor.jl")
export ChainMotor, dxdt_chainmotor

include("chain.jl")
include("disc_sim.jl")
export Chain, dxdt_chain, chain_odestep_DAE, chain_odestep_ODE

export dxdt_chain!, chainmotorDAE!, chainDAE!

include("misc.jl")
export getendpoint

end # module
