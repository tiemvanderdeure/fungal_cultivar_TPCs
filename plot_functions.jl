using CairoMakie, ColorSchemes
# convenience function
function quants_from_rates(f, ts, var; quantiles = [0.025, 0.5, 0.975])
    rate_est = f.(ts, eachrow(var)')
    reduce(hcat, quantile.(eachrow(rate_est), Ref(quantiles)))
end

function plot_species_mean!(ax, ts, tempmin, tempmax, chain, f; color = :black, linewidth = 3, transparency = 0.2, quantiles = [0.025, 0.5, 0.975], kw...)
    var_sp = Array(group(chain, :var_sp))
    sp_quants = quants_from_rates(f.defaults.f, ts, var_sp; quantiles = quantiles)
    ts_C = ts .- 273.15
    band!(ax, ts_C[tempmin:tempmax], sp_quants[1, tempmin:tempmax], sp_quants[3, tempmin:tempmax]; color = (color, transparency))
    lines!(ax, ts_C[101:301], sp_quants[2, 101:301]; color, linewidth, kw...)
    lines!(ax, ts_C[1:tempmin], sp_quants[2, 1:tempmin]; color, linewidth, linestyle = :dash)
    lines!(ax, ts_C[301:end], sp_quants[2, 301:end]; color, linewidth, linestyle = :dash)

    return maximum(sp_quants[3, :])
end

# main function
function plot_chain!(
    ax, metad, f, mychn, gq; 
    n_params, cols = ColorSchemes.Dark2_4,
    quantiles = [0.025, 0.5, 0.975]
    )
    ts_C = 0:0.1:50
    ts = ts_C .+ 273.15

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
    sp_quants = map(1:ncol) do i
        quants_from_rates(f.defaults.f, ts, var_co[:, (i*n_params - (n_params - 1)): i*n_params])
    end
    # lowest/highest temperature where there is any growth
    tempmin = max(min(minimum(map(s -> findfirst(s[2,:] .> 0.01), sp_quants)) - 5, 80), 1)
    tempmax = min(max(maximum(map(s -> findlast(s[2,:] .> 0.01), sp_quants)) + 5, 320), 501)

    for i in 1:ncol
        col = cols[i]
        lines!(ax, ts_C[101:301], sp_quants[i][2, 101:301], color = col)
        lines!(ax, ts_C[1:101], sp_quants[i][2, 1:101], color = col, linestyle = :dash)
        lines!(ax, ts_C[301:end], sp_quants[i][2, 301:end], color = col, linestyle = :dash)
    end

    xlims!(ax, ts_C[tempmin], ts_C[tempmax])

    # plot the central estimate
    upperlimit = plot_species_mean!(ax, ts, tempmin, tempmax, mychn, f)
    
    # Plot individual replicates as points
    xoffsets = rand(-0.5:0.001:0.5, n_indiv)
    Makie.scatter!(ax, metad[:, :temperature] .+ xoffsets, qrate_indiv[2, :], color = cols[metad[:, :colony_id]])
    Makie.errorbars!(
        ax, metad[:, :temperature] .+ xoffsets,
        qrate_indiv[2, :] , qrate_indiv[2, :] .- qrate_indiv[1, :], qrate_indiv[3, :] .- qrate_indiv[2, :],
        color = cols[metad[:, :colony_id]]
    )
    return upperlimit
end
