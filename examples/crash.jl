# This makes some versions of julia crash
using ChainDyn
chain = ChainMotor(Float64)
uconst = 10.0
foo(x) = dxdt_chain!(x, uconst, chain)
foo(randn(43))
