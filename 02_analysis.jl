include("bayesian_model.jl")
include("data.jl")
using Turing, ArviZ, Serialization, ThreadsX, Statistics, Tables

###### Reconstructing and selecting chains #######
# read in all the chains, if not loaded in the REPL already
#=
chains = map(models) do model
    map(species) do sp
        open("chains/$(model.shortname)/$(sp).jls", "r") do f
            deserialize(f)
        end
    end
end
=#
# Throw out all chains that converged
chains_converged = map(chains) do chns
    chns_converged = map(chns) do chn
        converged_chains_indices = map(1:5) do i
            rhat_deviation(chn[:, :, i]) |> mean < 0.03 # select chains with good convergence
        end
        chn[:, :, converged_chains_indices]
    end
end

# How many chains remaining?
chn_size = map(chains_converged) do chns
    map(chns) do chn
        size(chn)[3]
    end
end

# model functions to get generated quantities & assess model fit (WAIC)
fs = map(models) do model
    map(ds, metads) do d, metad
        func_one_species_indep(
            d.time, d.length, temperatures_K, metad.temperature_id, metad.colony_id, d.id;
            f = model.f, prior_mus = model.prior_mus, prior_std = model.prior_std, dist = MvNormal
        )
    end
end

# generated quantities
gqs = map((f,c) -> generated_quantities.(f, c), fs, chains_converged)

########## Summary statistics ############
# Save summary statistics
summaries = map(models, fs, chains_converged) do model, f_indep, chns
    summaries = ThreadsX.map(species, f_indep, chns) do sp, f, chn
        n_param = length(f.defaults.prior_mus)
        n_col = maximum(f.args.col_id)
        n_indiv = maximum(f.args.id)

        var_sp = group(chn, :var_sp)
        params_to_save = [var_sp.name_map.parameters; [:sigma, :indiv_std, :log_density]]

        # loop over each param, extract, save as a NamedTuple and merge into one big NamedTuple
        summary = mapreduce(merge, params_to_save) do p
            slice = chn[p]
            name = Tuple(Symbol.("$(p)_$t" for t in ["mean", "2.5%", "97.5%"]))
            NamedTuple{name}([mean(slice); quantile(slice, [0.025, 0.975])])
        end

        # get log likelihood values to calculate WAIC
        log_likelihood_raw = Turing.pointwise_loglikelihoods(
            f, MCMCChains.get_sections(chn, :parameters)
        )
        log_likelihood = (; obs = cat(values(log_likelihood_raw)..., dims = 3))

        idata = from_mcmcchains(
            chn;
            coords=(; param = 1:n_param, col=1:n_col, indiv = 1:n_indiv),
            dims=(col_std = (:param,), var_sp = (:param,), var_co = (:param,:col), rate_indiv = (:indiv,), intercept = (:indiv,)),
            log_likelihood,
            library="Turing"
        )

        WAIC = waic(idata)

        summary = merge(summary, WAIC.estimates)

        return merge((; species = sp), summary)
    end
    #CSV.write("tables/$(model.shortname)_summary_stats.csv", summaries)
    summaries
end    

summaries_all = map(summaries, models) do sum, model
    DataFrame(map(s -> merge((; model = model.shortname), s), sum))
end

summaries_all = vcat(summaries_all..., cols = :union)

CSV.write("tables/summaries_all.csv", summaries_all)

# Save relevant parts of the chains
map(models, outputs_converged) do model, (fs, chns, gqs)
    map((s, c) -> CSV.write("chains/$(model.shortname)/$s/params.csv", group(c, :var_sp)), species, chns)
    map((s, c) -> CSV.write("chains/$(model.shortname)/$(s)/colonies.csv", group(c, :var_co)), species, chns)
    map((s, c) -> CSV.write("chains/$(model.shortname)/$(s)/goodness_of_fit.csv", c[[:sigma, :indiv_std, :log_density]]), species, chns)
end    

# inverse structure - one csv for each specues with all models
summaries_by_model = map(enumerate(species)) do (i, s)
    summary_by_model = [sum[i] for sum in summaries]
    summary_by_model = map(summary_by_model, models) do s, model
        Base.structdiff(merge((; model = model.shortname), s), NamedTuple{(:species, )})
    end
    CSV.write("tables/$(s)_summary_stats.csv", summary_by_model)
    summary_by_model
end

###### Manuscript tables ########
# Table S2
gaussian_chns = chains_converged[1]
mod_gaussian_chns = chains_converged[2]
weibull_chns = chains_converged[3]
deutsch_chns = chains_converged[4]

function mat_to_quantile(mat)
    low, high = round.(quantile(mat, [0.025, 0.975]); sigdigits = 3)
    m = round(mean(mat); sigdigits = 3)
    return "$m (95%CI $low-$high)"
end
mat_to_quantile(::Nothing) = ""

gaussian_results = map(gaussian_chns) do chn
    topt = chn[Symbol("var_sp[1]")] .- 273.15
    rmax = exp.(chn[Symbol("var_sp[2]")]) .* 10
    std = exp.(chn[Symbol("var_sp[3]")])
    CTmin = topt .- 3 .* std
    CTmax = topt .+ 3 .* std
    breadth = CTmax .- CTmin
    post = (; CTmin, topt, CTmax, rmax, breadth)
    map(mat_to_quantile, post)
end |> DataFrame

mod_gaussian_results = map(mod_gaussian_chns) do chn
    topt = chn[Symbol("var_sp[1]")] .- 273.15
    rmax = exp.(chn[Symbol("var_sp[2]")]) .* 10
    params = (; CTmin = nothing, topt, CTmax = nothing, rmax, breadth = nothing)
    map(mat_to_quantile, params)
end |> DataFrame

weibull_results = map(weibull_chns) do chn
    topt = chn[Symbol("var_sp[1]")] .- 273.15
    rmax = exp.(chn[Symbol("var_sp[2]")]) .* 10
    params = (; CTmin = nothing, topt, CTmax = nothing, rmax, breadth = nothing)
    map(mat_to_quantile, params)
end |> DataFrame

deutsch_results = map(deutsch_chns) do chn
    topt = chn[Symbol("var_sp[1]")] .- 273.15
    CTmax = chn[Symbol("var_sp[2]")] .- 273.15
    rmax = exp.(chn[Symbol("var_sp[3]")]) .* 10
    std = exp.(chn[Symbol("var_sp[4]")])
    CTmin = topt .- 3 .* std
    breadth = CTmax .- CTmin
    params = (; CTmin, topt, CTmax, rmax, breadth)
    map(mat_to_quantile, params)
end |> DataFrame

table_S2 = vcat(gaussian_results, mod_gaussian_results, weibull_results, deutsch_results)

table_S2.model = repeat(getindex.(models, :shortname), inner = 6)
table_S2.ant_species = repeat(replace.(species, "_" => " "), outer = 4)
table_S2.farming_strategy = ifelse.(in.(table_S2.ant_species, Ref(["Cyphomyrmex costatus", "Apterostigma dentigerum"])), "Above", "Below")
select!(table_S2, [:ant_species, :farming_strategy, :model, :CTmin, :topt, :CTmax, :rmax, :breadth]) # change the order
sort!(table_S2, [:farming_strategy, :ant_species, :model])
CSV.write("tables/table_S2.csv", table_S2)

# deutsch only, for table 1
table_1 = filter(x -> x.model == "Deutsch", table_S2)
select!(table_1, Not(:model))
CSV.write("tables/table_1.csv", table_1)

###### Pairwise analyses ######
deutsch_topt = getindex.(deutsch_chns, Symbol("var_sp[1]"))
deutsch_ctmax = getindex.(deutsch_chns, Symbol("var_sp[2]"))
deutsch_rmax = getindex.(deutsch_chns, Symbol("var_sp[3]"))
deutsch_std = getindex.(deutsch_chns, Symbol("var_sp[4]"))
deutsch_ctmin = map((t, s) -> t .- 3 .* exp.(s), deutsch_topt, deutsch_std)
deutsch_breadth = map((max, min) -> max .- min, deutsch_ctmax, deutsch_ctmin)

weibull_topt = getindex.(weibull_chns, Symbol("var_sp[1]"))
weibull_rmax = getindex.(weibull_chns, Symbol("var_sp[2]"))
gaussian_rmax = getindex.(gaussian_chns, Symbol("var_sp[2]"))
gaussian_std = getindex.(gaussian_chns, Symbol("var_sp[3]"))
gauusian_topt = getindex.(gaussian_chns, Symbol("var_sp[1]"))


all_data = [deutsch_ctmax, deutsch_rmax, deutsch_ctmin, deutsch_breadth]
data_names = ["deutsch_ctmax", "deutsch_rmax", "deutsch_ctmin", "deutsch_breadth"]
for (data, name) in zip(all_data, data_names)
    matrix = [mean(d[:,1:min(size(d,2), size(d2,2))] .< d2[:,1:min(size(d,2), size(d2,2))]) for d in data, d2 in data]
    CSV.write("tables/$name.csv", Tables.table(matrix; header = species))
end

###### Correlations ######
deutsch_topt_mat = hcat(vec.(getindex.(deutsch_topt, :, Ref(1:4)))...)
deutsch_rmax_mat = exp.(hcat(vec.(getindex.(deutsch_rmax, :, Ref(1:4)))...))
deutsch_breadth_mat = hcat(vec.(getindex.(deutsch_breadth, :, Ref(1:4)))...)

topt_rmax = cor.(eachrow(deutsch_topt_mat), eachrow(deutsch_rmax_mat))
breadth_topt = cor.(eachrow(deutsch_breadth_mat), eachrow(deutsch_topt_mat))
breadth_rmax = cor.(eachrow(deutsch_breadth_mat), eachrow(deutsch_rmax_mat))

topt_rmax_prob = mean(topt_rmax .> 0)
breadth_topt_prob = mean(breadth_topt .> 0)
breadth_rmax_prob = mean(breadth_rmax .> 0)

# calculate bayes factors - since we know that the prior probability is 0.5
topt_rmax_bf = (1-topt_rmax_prob)/topt_rmax_prob
breadth_topt_bf = (1-breadth_topt_prob)/breadth_topt_prob
breadth_rmax_bf = (1-breadth_rmax_prob)/breadth_rmax_prob

topt, rmax, breadth = first(zip(eachrow(deutsch_topt_mat), eachrow(deutsch_rmax_mat), eachrow(deutsch_breadth_mat)))

using GLM
lms = map(eachrow(deutsch_topt_mat), eachrow(deutsch_rmax_mat), eachrow(deutsch_breadth_mat)) do topt, rmax, breadth
    nt = (; topt, rmax, breadth)
    rmax_topt = lm(@formula(rmax ~ topt), nt)
    breadth_topt = lm(@formula(breadth ~ topt), nt)
    breadth_rmax = lm(@formula(rmax ~ breadth), nt)
    return (; rmax_topt, breadth_topt, breadth_rmax)
end |> Tables.columntable

# Numbers mentioned in-text
correlations = map((;topt_rmax, breadth_topt, breadth_rmax)) do c
    [mean(c); quantile(c, [0.025, 0.975])...]
end
CSV.write("tables/correlations.csv", correlations)

####### Plots ###########
## Plot
include("plot_functions.jl")

species_order = [3,1,6,5,4,2]

figs = map(models, chains_converged, gqs, fs) do model, chns, gq, f
    # Figure, axes, and title
    fig = Figure(size  = (700, 700))
    figure_axes = [
        Axis(
            fig[i,j], 
            ylabel = "Growth rate (cm/day)", xlabel = "Temperature (°C)", titlefont = :bold_italic,
            limits = (nothing, nothing, 0, nothing),
            xgridvisible = false, ygridvisible = false
        ) 
        for j in 1:2, i in 1:3
    ]
    Label(fig[0, :], model.shortname, fontsize = 25, font = :bold)

    # Plot chains
    map(species_order, 1:6) do i, j
        figure_axes[j].title = replace(species[i], "_" => " ")
        upperlimit = plot_chain!(figure_axes[j], metads[i], f[i], chns[i], gq[i], n_params = f[i].defaults.n_params)
    end

    save("plots/$(model.shortname).pdf", fig; px_per_unit = 3)
    return fig
end

save("plots/figure_03.pdf", figs[4])

### Figure with mean estimates for A. dentigerum for all models
ts_c = 0.0:0.1:50.0
ts = ts_c .+ 273.15
fig = Figure(size = (800, 600))
ax = Axis(fig[1,1], limits = (minimum(ts) .- 273.15, maximum(ts) .- 273.15, 0., 0.12))

s_id = 2
for i in eachindex(models)
    chn = chains_converged[i][s_id]
    f = fs[i][s_id]
    gq = gqs[i][s_id]
    plot_species_mean!(ax, ts, 50, 450, chn, f; color = ColorSchemes.Dark2_7[i], linewidth = 3, transparency = 0, label = models[i].shortname)
end

axislegend(ax, position = :lt)

# Plot mean estimates for growth rate for each replicate
n_indiv = size(metads[s_id])[1]
gq = vec(gqs[1][s_id])
mean_rate_indiv = map(1:n_indiv) do i
    mean(map(g -> g.rate_indiv[i], gq))
end
scatter!(ax, metads[s_id][:, :temperature], mean_rate_indiv, color = :black)
fig

save("plots/figure_S1.pdf", fig)


### Figure with all species for the Deutsch model
fig = Figure(size = (800, 600))
ax = Axis(fig[1,1], limits = (8, 43, 0., 0.12))
chns_deutsch = chains_converged[4]
f_deutsch = fs[4]
for i in 1:6
    chn = chns_deutsch[i]
    f = f_deutsch[i]
    plot_species_mean!(
        ax, ts, 50, 450, chn, f; 
        color = ColorSchemes.Dark2_7[i], linewidth = 4, transparency = 0, label = replace(species[i], "_" => " "))
end
axislegend(ax, position = :lt)
save("plots/figure_04.pdf", fig)

gqs_converged = map((f,c) -> generated_quantities(f, c), fs, chns_converged)
return (fs, chns_converged, gqs_converged)

## Correlation plots
deutsch_topt_C = map(t -> t .- 273.15, deutsch_topt)
deutsch_rmax_exp = map(t -> exp.(t), deutsch_rmax)
is_aboveground = in.(species, Ref(["Cyphomyrmex_costatus", "Apterostigma_dentigerum"]))

posteriors = (topt = deutsch_topt_C, rmax = deutsch_rmax_exp, breadth = deutsch_breadth)
posterior_means = map(posteriors) do p
    mean.(p)
end
posterior_uncertainty = map(posteriors,posterior_means) do p, means
    mat = hcat(quantile.(p, Ref([0.025, 0.975]))...)
    [means .- mat[1,:] mat[2,:] .- means]
end

topts = (;topt = range(extrema(deutsch_topt_mat)...; length = 100))
breadths =(;breadth = range(extrema(deutsch_breadth_mat)...; length = 100))
toptc_C = topts.topt .- 273.15
rmax_topt_cor_mat = reduce(hcat, predict.(lms.rmax_topt, Ref(topts)))
breadth_topt_cor_mat = reduce(hcat, predict.(lms.breadth_topt, Ref(topts)))
breadth_rmax_cor_mat = reduce(hcat, predict.(lms.breadth_rmax, Ref(breadths)))

rmax_topt_cor_mean = mean(rmax_topt_cor_mat; dims = 2)[:,1]
rmax_topt_cor_quant = reduce(hcat, quantile.(eachrow(rmax_topt_cor_mat), Ref([0.025,0.2,0.8,0.975])))

breadth_topt_cor_mean = mean(breadth_topt_cor_mat; dims = 2)[:,1]
breadth_topt_cor_quant = reduce(hcat, quantile.(eachrow(breadth_topt_cor_mat), Ref([0.025,0.2,0.8,0.975])))

breadth_rmax_cor_mean = mean(breadth_rmax_cor_mat; dims = 2)[:,1]
breadth_rmax_cor_quant = reduce(hcat, quantile.(eachrow(breadth_rmax_cor_mat), Ref([0.025,0.2,0.8,0.975])))

f = Figure(size = (1700, 500))
ax1 = Axis(f[1,1],
    xlabel = "Topt (°C)", ylabel = "rmax (cm/day)",
    xgridvisible = false, ygridvisible = false)
ax2 = Axis(f[1,2],
    xlabel = "Topt (°C)", ylabel = "Thermal tolerance breadth (°C)",
    xgridvisible = false, ygridvisible = false)
ax3 = Axis(f[1,3],
    xlabel = "Thermal tolerance breadth (°C)", ylabel = "rmax (cm/day)",
    xgridvisible = false, ygridvisible = false)

color = is_aboveground
marker = [:circle, :rect, :diamond, :utriangle, :dtriangle, :pentagon]
colormap = [:blue, :red]
markersize = 18
sc = scatter!(ax1, 
    posterior_means.topt, posterior_means.rmax;
    color, colormap, markersize, marker)
Makie.errorbars!(ax1, 
    posterior_means.topt, posterior_means.rmax,
    posterior_uncertainty.topt[:,1], posterior_uncertainty.topt[:,2];
    direction = :x, color, colormap)
Makie.errorbars!(ax1, 
    posterior_means.topt, posterior_means.rmax,
    posterior_uncertainty.rmax[:,1], posterior_uncertainty.rmax[:,2];
    color, colormap)


scatter!(ax2, 
    posterior_means.topt, posterior_means.breadth;
    color, colormap, markersize, marker)
Makie.errorbars!(ax2, 
    posterior_means.topt, posterior_means.breadth,
    posterior_uncertainty.topt[:,1], posterior_uncertainty.topt[:,2];
    direction = :x, color, colormap)
Makie.errorbars!(ax2, 
    posterior_means.topt, posterior_means.breadth,
    posterior_uncertainty.breadth[:,1], posterior_uncertainty.breadth[:,2];
    color, colormap)
   
scatter!(ax3, 
    posterior_means.breadth, posterior_means.rmax;
    color, colormap, markersize, marker)
Makie.errorbars!(ax3, 
    posterior_means.breadth, posterior_means.rmax,
    posterior_uncertainty.breadth[:,1], posterior_uncertainty.breadth[:,2];
    direction = :x, color, colormap)
Makie.errorbars!(ax3, 
    posterior_means.breadth, posterior_means.rmax,
    posterior_uncertainty.rmax[:,1], posterior_uncertainty.rmax[:,2];
    color, colormap)

legend_aboveground = [MarkerElement(; marker = marker[i], color = colormap[2], strokecolor = :transparent, markersize) for i in findall(is_aboveground)]
legend_belowground = [MarkerElement(; marker = marker[i], color = colormap[1], strokecolor = :transparent, markersize) for i in findall(.!is_aboveground)]
species_label = replace.(species, "_" => " ")
Legend(f[1,4], [legend_aboveground, legend_belowground], [species_label[is_aboveground], species_label[.!is_aboveground]], ["Aboveground", "Belowground"])
f # need to call here to generate ax1.finalllimit[]

for (i, (ax, x, y_m, y_q)) in enumerate(zip(
    (ax1, ax2, ax3), (toptc_C, toptc_C, breadths.breadth), 
    (rmax_topt_cor_mean, breadth_topt_cor_mean, breadth_rmax_cor_mean), 
    (rmax_topt_cor_quant, breadth_topt_cor_quant, breadth_rmax_cor_quant)))
    if i > 1
        limits!(ax, ax.finallimits[])
        lines!(ax, x, y_m; color = :black)
        lines!(ax, x, y_q[1,:]; color = :black, linestyle = :dot)
        lines!(ax, x, y_q[2,:]; color = :black, linestyle = :dash)
        lines!(ax, x, y_q[3,:]; color = :black, linestyle = :dash)
        lines!(ax, x, y_q[4,:]; color = :black, linestyle = :dot)
    end
end

save(f, "plots/figure_05.pdf")