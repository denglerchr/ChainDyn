struct Foo{T}
    M::Array{T, 2}
end

function dostuff!(x, u::T, foo::Foo{T}) where {T}
    Msym = Symmetric(foo.M, :U)
    return 1.0
end

# Run example
myfoo = Foo(Array{Float64}(I, 3, 3))
myu = 1.0

test(x) = dostuff!(x, myu, myfoo)

# "dostuff!(1, myu, myfoo)" works for me
test(1)
