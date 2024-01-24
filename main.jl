include("data.jl")
include("bayesian_model.jl"); include("plot_functions.jl")
include("response_functions.jl")

Turing.setadbackend(:forwarddiff)
Turing.setrdcache(true)
using Random

# Select models and priors
prior_topt = 25 + 273.15
models = [
    (f = gaussian_from_log_v, shortname = "Gaussian", prior_mus = [prior_topt, -2.5, 2.], prior_std = [5, 1., 1.], adbackend = :reversediff),
    (f = mod_gaussian_from_log_v, shortname = "Modified Gaussian", prior_mus = [prior_topt, -2.5, 2., 0.7], prior_std = [5, 1., 1., 0.5], adbackend = :reversediff),
    (f = weibull_from_logv, shortname = "Weibull", prior_mus = [prior_topt, -2.5, 4, -1], prior_std = [5, 1., 2., 2.], adbackend = :forwarddiff),
    (f = pawar_from_logv, shortname = "Pawar", prior_mus = [prior_topt, -3.0, 0.5, 1.5], prior_std = [5, 0.5, 1., 1.], adbackend = :forwarddiff),
    (f = deutsch_from_logv, shortname = "Deutsch", prior_mus = [prior_topt, prior_topt+10, -2.5, 1.], prior_std = [5., 5., 0.5, 1.], adbackend = :forwarddiff),
]

# Run once, make directories
for model in models 
    dir = "chains/$(model.shortname)"
    if !isdir(dir)
        mkdir(dir)
    end	
        for s in species
        dir = "chains/$(model.shortname)/$s"
        if !isdir(dir)
            mkdir(dir)
        end	
    end
end

# Prepare data
ds = map(1:6) do sp_id
    d = filter(x -> x.species_id == sp_id, dlong)
    d.colony_id = StatsBase.denserank(d.colony_id)
    d.id = StatsBase.denserank(d.id)
    return d
end

metads = map(d -> d[d.time .== 0, [:id, :temperature_id, :colony_id, :T_K, :temperature]], ds)

outputs = map(models) do model
    println("starting on model $(model.shortname)")
    Turing.setadbackend(model.adbackend)

    # Run the models
    fs = map(ds, metads) do d, metad
        func_one_species(
            d.time, d.length, temperatures_K, metad.temperature_id, metad.colony_id, d.id;
            f = model.f, prior_mus = model.prior_mus, prior_std = model.prior_std, dist = MvNormal
        )
    end

    chns = map(fs) do f
        Random.seed!(0) # for reproducibility
        sample(f, NUTS(), MCMCThreads(), 2000, 5)
    end
    #chns = map(f -> sample(f, NUTS(), 1500), fs)
    gqs = map((f,c) -> generated_quantities(f, c), fs, chns)

    return (fs, chns, gqs)
end

# A few chains went haywire. Throw those out and select all chains that did converge.
outputs_converged = map(outputs) do (fs, chns, gqs)
    chns_converged = map(chns) do chn
        converged_chains_indices = map(1:5) do i
            rhat_deviation(chn[:, :, i]) |> mean < 0.03 # select chains with good convergence
        end
        chn[:, :, converged_chains_indices]
    end

    gqs_converged = map((f,c) -> generated_quantities(f, c), fs, chns_converged)
    return (fs, chns_converged, gqs_converged)
end

chn_size = map(outputs_converged) do (fs, chns, gqs)
    map(chns) do chn
        size(chn)[3]
    end
end

# Save relevant parts of the chains
map(models, outputs_converged) do model, (fs, chns, gqs)
    map((s, c) -> CSV.write("chains/$(model.shortname)/$s/params.csv", group(c, :var_sp)), species, chns)
    map((s, c) -> CSV.write("chains/$(model.shortname)/$(s)/colonies.csv", group(c, :var_co)), species, chns)
    map((s, c) -> CSV.write("chains/$(model.shortname)/$(s)/goodness_of_fit.csv", c[[:sigma, :indiv_std, :log_density]]), species, chns)
end    

## Plot
    map(models, outputs_converged) do model, (fs, chns, gqs)
    # Figure, axes, and title
    fig = Figure(size  = (700, 700))
    figure_axes = [Axis(fig[i,j], ylabel = "growth (cm/day)", xlabel = "temperature (°C)", titlefont = :bold_italic) for j in 1:2, i in 1:3]
    Label(fig[0, :], model.shortname, fontsize = 25, font = :bold)

    # Plot chains
    for i in 1:6
        figure_axes[i].title = replace(species[i], "_" => " ")
        plot_chain!(figure_axes[i], metads[i], fs[i], chns[i], gqs[i], n_params = fs[i].defaults.n_params)
    end

    save("plots/$(model.shortname).png", fig; px_per_unit = 3)
end


## Plot
map(models, outputs_converged) do model, (fs, chns, gqs)
    summary_stats = map(chns, species) do chn, species # oop over each species and combine
        var_sp = group(chn, :var_sp)
        params_to_save = [var_sp.name_map.parameters; [:sigma, :indiv_std, :log_density]]

        # loop over each param, extract, save as a NamedTuple and merge into one big NamedTuple
        params_summary = mapreduce(merge, params_to_save) do p
            slice = chn[p]
            name = Tuple(Symbol.("$(p)_$t" for t in ["mean", "2.5%", "97.5%"]))
            NamedTuple{name}([mean(slice); quantile(slice, [0.025, 0.975])])
        end

        return merge((; species = species), params_summary)

    end
    CSV.write("tables/$(model.shortname)_summary_stats.csv", summary_stats)
end


using ColorSchemes
colscheme = ColorSchemes.Dark2_7

ts = 280:0.1:310
fig = Figure(size = (800, 600))
ax = Axis(fig[1,1], limits = (minimum(ts) .- 273.15, maximum(ts) .- 273.15, 0., 0.12))

s_id = 2
for i in 1:5
    (fs, chns, gqs) = outputs_converged[i]
    plot_species_mean!(ax, ts, chns[s_id], fs[s_id]; color = colscheme[i], linewidth = 3, transparency = 0, label = models[i].shortname)
end

axislegend(ax, position = :lt)

# Plot mean estimates for growth rate for each replicate
n_indiv = size(metads[s_id])[1]
gq = vec(outputs_converged[1][3][s_id])
mean_rate_indiv = map(1:n_indiv) do i
    mean(map(g -> g.rate_indiv[i], gq))
end
scatter!(ax, metads[s_id][:, :temperature], mean_rate_indiv, color = :black)
fig

save("plots/all_models_A_dentigerum.png", fig)
