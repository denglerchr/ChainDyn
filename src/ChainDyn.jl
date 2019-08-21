module ChainDyn

using LinearAlgebra

abstract type AbstractChain end # general chain type

include("forcevec20.jl")
include("massmat20.jl")
include("chainmotor.jl")
export ChainMotor, dxdt_chainmotor

include("chain.jl")
include("disc_sim.jl")
export Chain, dxdtchain, chain_odestep

export dxdt_chain!, chainmotorDAE!, chainDAE!

end # module
