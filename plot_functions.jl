using CairoMakie
# convenience function
function quants_from_rates(f, ts, var; quantiles = [0.025, 0.5, 0.975])
    rate_est = f.(ts, eachrow(var)')
    reduce(hcat, quantile.(eachrow(rate_est), Ref(quantiles)))
end

function plot_species_mean!(ax, ts, chain, f; color = :black, linewidth = 3, transparency = 0.2, quantiles = [0.025, 0.5, 0.975], kw...)
    var_sp = Array(group(chain, :var_sp))
    sp_quants = quants_from_rates(f.defaults.f, ts, var_sp; quantiles = quantiles)
    ts_C = ts .- 273.15
    lines!(ax, ts_C, sp_quants[2, :]; color = color, linewidth = linewidth, kw...)
    band!(ax, ts_C, sp_quants[1, :], sp_quants[3, :]; color = (color, transparency))
end

# main function
function plot_chain!(
    ax, metad, f, mychn, gq; 
    n_params, ts = collect(281:0.1:305), cols = [:red, :blue, :green, :orange],
    quantiles = [0.025, 0.5, 0.975]
    )
    ts_C = ts .- 273.15

    n_indiv = size(metad)[1]
    var_co = Array(group(mychn, :var_co))
    rate_indiv = map(1:n_indiv) do i
        map(vec(gq)) do g
            g.rate_indiv[i]
        end
    end

    qrate_indiv = mapreduce(x -> quantile(x, quantiles), hcat, rate_indiv)

    ncol = Int(size(var_co)[2] / n_params)
    # plot each colony
    for i in 1:ncol
        col = cols[i]
        sp_quants = quants_from_rates(f.defaults.f, ts, var_co[:, (i*n_params - (n_params - 1)): i*n_params])
        
        lines!(ax, ts_C, sp_quants[2, :], color = col)
        #band!(ax, ts, sp_quants[1, :], sp_quants[3, :], color = (col, 0.2))
    end
    
    # plot the central estimate
    plot_species_mean!(ax, ts, mychn, f)
    
    xoffsets = rand(-0.5:0.001:0.5, n_indiv)
    Makie.scatter!(ax, metad[:, :temperature] .+ xoffsets, qrate_indiv[2, :], color = cols[metad[:, :colony_id]])
    Makie.errorbars!(
        ax, metad[:, :temperature] .+ xoffsets,
        qrate_indiv[2, :] , qrate_indiv[2, :] .- qrate_indiv[1, :], qrate_indiv[3, :] .- qrate_indiv[2, :],
        color = cols[metad[:, :colony_id]]
    )
end
