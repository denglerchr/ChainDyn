# Coulomb+viscous friction
function carfriction(v::T, Fe::T) where {T}
    Fs = T(76.697)
    cv = T(92.97)
    if  abs(v)>=eps()
        return T(sign(v)*(Fs+cv*abs(v)))
    elseif v<eps() && abs(Fe)<Fs
        return T(Fe)
    else
        return T( sign(Fe)*Fs )
    end
end

# Differentiable approximation of the above
function carfriction_approx(v::T, Fe::T) where {T}
    Fs = T(76.697)
    cv = T(92.97)
    Fcoulomb = Fs*tanh(T(20)*v)
    Fviscous = cv*v
    return Fcoulomb + Fviscous
end
