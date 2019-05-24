module ChainDyn

using LinearAlgebra

include("forcevec20.jl")
include("massmat20.jl")
include("chainmotor.jl")
export ChainMotor, dxdt_chainmotor

include("chain.jl")
export Chain, dxdtchain

export dxdt_chain!, chainmotorDAE!, chainDAE!

end # module
