using CSV, DataFrames, StatsBase
d = CSV.read("data.csv", DataFrame)
times = collect(0:7:28)
d.id = 1:size(d,1)

dlong = stack(d, filter(x -> occursin("area", x), names(d)); variable_name = :time, value_name = :area)
dlong.time = parse.(Int, getindex.(split.(dlong.time, "_d"), 2)) # time is an integer

dropmissing!(dlong, :area) # drop missing area values (infected)

dlong.length = sqrt.(dlong.area ./ pi) # get radius
dlong.species_id = StatsBase.denserank(dlong.species)
dlong.colony_id = StatsBase.denserank(dlong.colony)
dlong.temperature_id = StatsBase.denserank(dlong.temperature)
dlong.T_K = dlong.temperature .+ 273.15

temp_temp_id = unique(dlong[!, [:temperature, :temperature_id]])
col_id_species_id = Vector(sort(unique(dlong[!, [:colony_id, :species_id]]), :colony_id).species_id)

# some metadata
species = [dlong[findfirst(dlong.species_id .== i), :species] for i in 1:6]
temperatures_K = sort(unique(dlong.T_K))

# Prepare data
ds = map(1:6) do sp_id
    d = filter(x -> x.species_id == sp_id, dlong)
    d.colony_id = StatsBase.denserank(d.colony_id)
    d.id = StatsBase.denserank(d.id)
    return d
end

metads = map(d -> d[d.time .== 0, [:id, :temperature_id, :colony_id, :T_K, :temperature]], ds)
