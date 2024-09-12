include("response_functions.jl")
## Import some packages
using Turing, Distributions, StatsBase, ReverseDiff

@model function func_one_species(
    time, obs, temperatures, temp_id, col_id, id,
    n = maximum(id), ncol = maximum(col_id);
    f, prior_mus, prior_std, n_params = size(prior_mus)[1], dist = MvNormal
)

    # hyperparameters
    sigma ~ Exponential(1)
    indiv_std ~ Exponential(1) #prior on variation between individuals within a colony
    col_std ~ filldist(Exponential(1), n_params) #prior on varation between colonies

    # Species level
    var_sp ~ dist(prior_mus, prior_std)

    # Colony level
    var_co ~ filldist(dist(var_sp, col_std), ncol)

    Turing.@addlogprob! (added_logpdf(f, var_sp)) # to prevent invalid species means

    # rates for each temperature and each colony
    est_replicates = f.(temperatures, eachcol(var_co)')
    log_est_replicates = NaNMath.log.(est_replicates)

    # individual level
    est_indiv = map((t, c) -> getindex(log_est_replicates, t, c), temp_id, col_id)
    rate_indiv ~ MvLogNormal(est_indiv, indiv_std)
    intercept ~ filldist(Normal(0.2725, 0.2), n) # size at time 0

    # sampling
    length_estimate = (rate_indiv[id] .* time .+ intercept[id])
    obs ~ MvNormal(length_estimate, sigma)

    return (; rate_indiv, est_replicates)
end

# The exact same model as above, but with sampling in a loop instead of with MvNormal.
# This is much slower, but observations appear independent (which they are), which is needed for estiamting model performance
@model function func_one_species_indep(
    time, obs, temperatures, temp_id, col_id, id,
    n = maximum(id), ncol = maximum(col_id);
    f, prior_mus, prior_std, n_params = size(prior_mus)[1], dist = MvNormal
)

    # hyperparameters
    sigma ~ Exponential(1)
    indiv_std ~ Exponential(1) #prior on variation between individuals within a colony
    col_std ~ filldist(Exponential(1), n_params) #prior on varation between colonies

    # Species level
    var_sp ~ dist(prior_mus, prior_std)

    # Colony level
    var_co ~ filldist(dist(var_sp, col_std), ncol)

    Turing.@addlogprob! (added_logpdf(f, var_sp)) # to prevent invalid species means

    # rates for each temperature and each colony
    est_replicates = f.(temperatures, eachcol(var_co)')
    log_est_replicates = NaNMath.log.(est_replicates)

    # individual level
    est_indiv = map((t, c) -> getindex(log_est_replicates, t, c), temp_id, col_id)
    rate_indiv ~ MvLogNormal(est_indiv, indiv_std)
    intercept ~ filldist(Normal(0.2725, 0.2), n) # size at time 0

    # sampling
    length_estimate = (rate_indiv[id] .* time .+ intercept[id])
    obs .~ Normal.(length_estimate, sigma)
    #for i in 1:n
    #    bv = id .== i
    #    obs[bv] ~ MvNormal(length_estimate[bv], sigma)
    #end

    return (; rate_indiv, est_replicates)
end

# Select models and priors
prior_topt = 25 + 273.15
models = [
    (f = gaussian_from_log_v, shortname = "Gaussian", prior_mus = [prior_topt, -2.5, 2.], prior_std = [5, 1., 1.], adbackend = AutoReverseDiff(true)),
    (f = mod_gaussian_from_log_v, shortname = "Modified Gaussian", prior_mus = [prior_topt, -2.5, 2., 0.7], prior_std = [5, 1., 1., 0.5], adbackend = AutoReverseDiff(true)),
    (f = weibull_from_logv, shortname = "Weibull", prior_mus = [prior_topt, -2.5, 4, -1], prior_std = [5, 1., 2., 2.], adbackend = AutoReverseDiff(false)),
    (f = deutsch_from_logv, shortname = "Deutsch", prior_mus = [prior_topt, prior_topt+7, -2.5, 1.], prior_std = [5., 3., 0.5, 1.], adbackend = AutoReverseDiff(false)),
]

# some tools
rhat_deviation(chain) = abs.(1 .- rhat(chain).nt[2])

# Some convergence check, to be improved
function converged(chain)
    absrh = rhat_deviation(chain)
    if mean(absrh) > 0.02 return false end
    if maximum(absrh) > 0.1 return false end
    return true
end
