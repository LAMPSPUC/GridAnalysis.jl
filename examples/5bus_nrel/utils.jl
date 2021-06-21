"""
    build_5_bus_matpower_DA(data_dir::AbstractString, case_file::AbstractString)

Builds base system for the 5bus NREL case (a.k.a NESTA case) from:
 - A matpower file containing grid information (case_file);
 - A file describing forecasts locations and details (forecasts_pointers_file);
 - A simple definition of reserves.
"""
function build_5_bus_matpower_DA(
    data_dir::AbstractString;
    case_file::AbstractString="case5_re_uc.m",
    FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
    forecasts_pointers_file::AbstractString=joinpath(
        FORECASTS_DIR, "timeseries_pointers_da_7day.json"
    ),
    add_reserves::Bool=true,
)
    case_file_path = joinpath(data_dir, case_file)
    sys = System(case_file_path)
    if add_reserves
        reserves = [
            VariableReserve{ReserveUp}("REG1", true, 5.0, 0.2),
            VariableReserve{ReserveUp}("REG2", true, 5.0, 0.4),
            VariableReserve{ReserveUp}("REG3", true, 5.0, 0.3),
            VariableReserve{ReserveUp}("REG4", true, 5.0, 0.4),
            VariableReserve{ReserveDown}("REG1d", true, 5.0, 0.2),
            VariableReserve{ReserveDown}("REG2d", true, 5.0, 0.4),
            VariableReserve{ReserveDown}("REG3d", true, 5.0, 0.3),
            VariableReserve{ReserveDown}("REG4d", true, 5.0, 0.4),
        ]
        contributing_devices = get_components(Generator, sys)
        for r in reserves
            add_service!(sys, r, contributing_devices)
        end
    end

    add_time_series!(sys, forecasts_pointers_file)

    return sys
end


"""
    build_5_bus_matpower_RT(; kwargs...)

Builds base system for the 5bus NREL case (a.k.a NESTA case) from:
 - A matpower file containing grid information (case_file);
 - A file describing forecasts locations and details (forecasts_pointers_file);
"""
function build_5_bus_matpower_RT(data_dir::AbstractString;
    case_file::AbstractString="case5_re_uc.m",
    FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
    forecasts_pointers_file::AbstractString=joinpath(
        FORECASTS_DIR, "timeseries_pointers_rt_7day.json"
    )
)
    case_file_path = joinpath(data_dir, case_file)
    sys = System(case_file_path)

    add_time_series!(sys, forecasts_pointers_file)
    transform_single_time_series!(sys, 1, Minute(5))

    return sys
end


"""
    prep_systems_UCED(system::System)

Duplicates the system to represent UC and ED for DA, transforming the time series
to the appropriate interval and horizon.
"""
function prep_systems_UCED(
    system::System;
    horizon_uc::Int=24,
    horizon_ed::Int=1,
    interval_uc::TimePeriod=Hour(24),
    interval_ed::TimePeriod=Hour(1),
)
    system_uc = system
    system_ed = deepcopy(system)

    transform_single_time_series!(system_uc, horizon_uc, interval_uc)
    transform_single_time_series!(system_ed, horizon_ed, interval_ed)

    return system_uc, system_ed
end

"""
    plot_price_curves(
        lmps_df::Dict{Any, Any}, period::Vector{Int64}, 
        bus_name::AbstractArray=["bus5"]
    )

Function to plot the price curve for the virtual offer bids. 
The 'bus_names' and 'periods' controls which buses and periods we want to include
in the plot, respectively.
"""

function plot_price_curves(lmps_df::Dict{Any, Any}, period::Vector{Int64}, bus_name::AbstractArray=["bus5"], node:: String="bus5")

    lmps_df = sort(lmps_df)
    indices=[]
    data=Array{Any}(nothing, (length(period),length(bus_name)+1,length(lmps_df))) #length(period)-size(lmps_df[0])[1]
    for (i, v) in enumerate(keys(lmps_df))
        for t =1:length(period) 
            data[t,1,i]=lmps_df[v][!,"DateTime"][period[t]]
            c=2
            for j in bus_name
                data[t,c,i]=lmps_df[v][!,j][period[t]]
                c=c+1
            end
        end
        indices = vcat(indices, v)
    end
   
    c=1
    for b=1: length(bus_name) 
        for t=1:length(period)
            if c==1
                plot(indices, data[t,b+1,:],label="hora"*string(period[t]-1)*"- "*string(bus_name[b]), legend=:outertopright)
            else
                plot!(indices, data[t,b+1,:],label="hora"*string(period[t]-1)*"- "*string(bus_name[b]), legend=:outertopright)
            end
            c=c+1
        end
    end

    plot!(title = "Price per Virtual Bid on "*node, ylabel = "Prices (\$/MWh)", xlabel = "Bid offers (p.u.)")
    #TODO: Change x axis, define bus_name of virtual bid 
end

"""
    plot_revenue_curves(
        lmps_df::Dict{Any, Any}, results_df::Dict{Any, Any}, market_simulator::MarketSimulator, period::Vector{Int64}, 
        bus_name::AbstractArray=["bus5"],
    )

Function to plot the revenue curve for the the virtual offer bids. 
The 'generator_name' defines which is the virtual generator that we want to plot it's results
and 'periods' controls which periods we want to include in the plot.
"""
function plot_revenue_curves(lmps_df::Dict{Any, Any}, results_df, market_simulator::MarketSimulator, period::Vector{Int64}, generator_name::String)
    
    lmps_df = sort(lmps_df)
    gen=get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_name=get_name(get_bus(gen))

    indices=[]
    data=Array{Any}(nothing, (length(period),2,length(lmps_df)))
    for (i, v) in enumerate(keys(lmps_df))
        for t =1:length(period)
            data[t,1,i]=lmps_df[v][!,"DateTime"][period[t]]
            variable_results = read_realized_variables(get_problem_results(results_df[v], "UC"), names=[:P__ThermalStandard])
            generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
            virtual_gen=generator_data[1][!,generator_name][[period[t]]][1] #name_generator
            data[t,2,i]=(lmps_df[v][[period[t]],bus_name][1])*virtual_gen #arrumar 
        end
        indices = vcat(indices, v)
    end

    c=1
    for t =1:length(period)
        if c==1
            plot(indices, data[t,2,:],label="hora"*string(period[t]-1), legend=:outertopright)
        else
            plot!(indices, data[t,2,:],label="hora"*string(period[t]-1), legend=:outertopright)
        end
        c=c+1
    end

    plot!(title = "Virtual Revenue per Offer on "*bus_name, ylabel = "Revenue (\$/MWh)", xlabel = "Bid offers (p.u)")
    #TODO: Change x axis 

end

"""
    plot_revenue_curves(
        lmps_df::Dict{Any, Any}, results_df::Dict{Any, Any}, market_simulator::MarketSimulator, period::Vector{Int64},
    )

Function to plot the virtual generation curve for the virtual offer bids. 
The 'generator_name' defines which is the virtual generator that we want to plot it's results
and 'periods' controls which periods we want to include in the plot.
"""

function plot_generation_curves(lmps_df::Dict{Any, Any}, results_df::Dict{Any, Any}, market_simulator::MarketSimulator, period::Vector{Int64}, generator_name::String)
    lmps_df=sort(lmps_df)
    gen=get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_name=get_name(get_bus(gen))

    indices=[]
    data=Array{Any}(nothing, (length(period),2,length(lmps_df)))
    for (i, v) in enumerate(keys(lmps_df))
        for t =1:length(period)
            data[t,1,i]=lmps_df[v][!,"DateTime"][period[t]]
            variable_results = read_realized_variables(get_problem_results(results_df[v], "UC"), names=[:P__ThermalStandard])
            generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
            virtual_gen=generator_data[1][!,generator_name][[period[t]]][1] 
            data[t,2,i]=virtual_gen #arrumar 
        end
        indices = vcat(indices, v)
    end

    c=1
    for t =1:length(period)
        if c==1
            plot(indices, data[t,2,:],label="hora"*string(period[t]-1), legend=:outertopright)
        else
            plot!(indices, data[t,2,:],label="hora"*string(period[t]-1), legend=:outertopright)
        end
        c=c+1
    end

    plot!(title = generator_name*" generation per Offer on "*bus_name, ylabel = "Generation(MWh)", xlabel = "Bid offers (p.u)")

end
