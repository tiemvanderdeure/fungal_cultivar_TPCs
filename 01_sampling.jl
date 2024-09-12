include("data.jl")
include("bayesian_model.jl")
using ThreadsX, Random, Serialization

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

chains = ThreadsX.map(models) do model 
    println("starting on model $(model.shortname)")
    # Run the models
    ThreadsX.map(zip(ds, metads, species)) do (d, metad, sp)
        # model function with data
        f = func_one_species(
            d.time, d.length, temperatures_K, metad.temperature_id, metad.colony_id, d.id;
            f = model.f, prior_mus = model.prior_mus, prior_std = model.prior_std, dist = MvNormal
        )
        
        # Sample the posterior
        chn = sample(Xoshiro(1234), f, NUTS(; adtype = model.adbackend), MCMCThreads(), 2000, 5)

        # Save the chain
        open("chains/$(model.shortname)/$(sp).jls", "w") do f
            serialize(f, chn)
        end
        return chn
    end
end

