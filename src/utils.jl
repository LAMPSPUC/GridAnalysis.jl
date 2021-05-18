"""
    fuel_type_mapping(system)

    return a mapping between the bus name and the fuel type for the given `system`.
"""
function fuel_type_mapping(system::System)
    generator_metadata =  [gen for gen in get_components(Generator, system)]

    fuel_type_map = Dict()
    for generator in generator_metadata 
        name = generator.name
        try 
            fuel_type_map[name] =  get_fuel(generator)
        catch
            # Assumes the bus name has format "<Fuel>Bus..."
            fuel_type_map[name] = first(split(name, "Bus"))
        end
    end

    return fuel_type_map
end
