function getendpoint(chain::AbstractChain, dtype = eltype(chain.q))
    # Constants
    l_link = dtype(1.178/20)

    # Calculate x coordicante of the endpoint
    Out_x::dtype = chain.q[1]
    Out_y::dtype = zero(dtype)
    for alpha in chain.q[2:end]
        Out_x += sin(alpha)*l_link
        Out_y -= cos(alpha)*l_link
    end
    return Out_x, Out_y
end

function getendvelocity(c::Chain, dtype::Type = eltype(c.u))
    l_link = dtype(1.178/20)
    v::dtype = c.u[1]
    for i = 2:length(c.u)
        v += c.u[i]*cos(c.q[i])*l_link
    end
    return v
end

function getendpoint(x::Vector{T}) where {T<:Number}
    @assert length(x) == 42
    c = Chain(T)
    c.q .= x[1:21]
    c.u .= x[22:end]
    return getendpoint(c, T)
end

function getendvelocity(x::Vector{T}) where {T<:Number}
    @assert length(x) == 42
    c = Chain(T)
    c.q .= x[1:21]
    c.u .= x[22:end]
    return getendvelocity(c, T)
end
