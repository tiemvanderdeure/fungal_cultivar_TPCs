using NaNMath

mod_gaussian(T, topt, a, b, c) = a * exp(-0.5*(abs(T - topt) / b)^c)
mod_gaussian_from_log_v(T, vars) = mod_gaussian(T, vars[1], exp(vars[2]), exp(vars[3]), exp(vars[4]))

gaussian(T, topt, a, b) = mod_gaussian(T, topt, a, b, 2)
gaussian_from_log_v(T, vars) = gaussian(T, vars[1], exp(vars[2]), exp(vars[3]))

function weibull(temp, topt, a, b, c)
    x = (temp - topt) / b + NaNMath.pow((c-1)/c, 1/c)
    if x > 0 
        o = a * NaNMath.pow((c-1)/c, (1-c)/c) * NaNMath.pow(x, c-1) * exp(-NaNMath.pow(x, c) + ((c-1)/c))
    else
        o = 0.
    end
    return o
end
weibull_from_logv(T, vars) = weibull(T, vars[1], exp(vars[2]), exp(vars[3]), exp(vars[3] + vars[4]))


function pawar(temp, topt, r_tref, e, eh; tref = 293.15)
    k = 8.62e-05
    r = r_tref * exp(-e/k * (1/temp - 1/tref)) / (1 + (e / (eh - e)) * exp(eh/k * (1/topt - 1/temp)))
    return r
end
pawar_from_logv(T, vars) = pawar(T, vars[1], exp(vars[2]), exp(vars[3]), exp(vars[4]))

# Add logpdf to prevent unvalid values for species mean
# by default 0
added_logpdf(f, x) = 0.

function deutsch(temp, topt, ctmax, rm, a)
    if temp > ctmax || topt > ctmax
        zero(temp)
    elseif temp > topt
        rm * (1 - ((temp - topt) / (topt - ctmax))^2)
    else
        gaussian(temp, topt, rm, a)
    end
end

deutsch_from_logv(T, vars) = deutsch(T, vars[1], vars[2], exp(vars[3]), exp(vars[4]))

# for pawar, eh > e
function added_logpdf(f::typeof(pawar_from_logv), x) 
    if x[4] < x[3] # if eh is lower than e, reject
        -Inf
    else
        0.
    end
end


# for deutsch, topt > ctmax
function added_logpdf(f::typeof(deutsch_from_logv), x) 
    if x[1] > x[2] # if topt is higher than ctmax, reject
        -Inf
    else
        0.
    end
end


# Fix what seems like an oversight in forwarddiff?
using ForwardDiff
using ForwardDiff: value, partials, Dual, isconstant, _mul_partials

ForwardDiff.@define_binary_dual_op(
    NaNMath.pow,
    begin
        vx, vy = value(x), value(y)
        expv = (NaNMath.pow)(vx, vy)
        powval = vy * (NaNMath.pow)(vx, vy - 1)
        if isconstant(y)
            logval = one(expv)
        elseif iszero(vx) && vy > 0
            logval = zero(vx)
        else
            logval = expv * NaNMath.log(vx)
        end
        new_partials = _mul_partials(partials(x), partials(y), powval, logval)
        return Dual{Txy}(expv, new_partials)
    end,
    begin
        v = value(x)
        expv = (NaNMath.pow)(v, y)
        if y == zero(y) || iszero(partials(x))
            new_partials = zero(partials(x))
        else
            new_partials = partials(x) * y * (NaNMath.pow)(v, y - 1)
        end
        return Dual{Tx}(expv, new_partials)
    end,
    begin
        v = value(y)
        expv = (NaNMath.pow)(x, v)
        deriv = (iszero(x) && v > 0) ? zero(expv) : expv*log(x)
        return Dual{Ty}(expv, deriv * partials(y))
    end
)
