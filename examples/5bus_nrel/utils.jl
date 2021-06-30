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
function build_5_bus_matpower_RT(
    data_dir::AbstractString;
    case_file::AbstractString="case5_re_uc.m",
    FORECASTS_DIR::AbstractString=joinpath(data_dir, "forecasts"),
    forecasts_pointers_file::AbstractString=joinpath(
        FORECASTS_DIR, "timeseries_pointers_rt_7day.json"
    ),
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

function plot_price_curves(
    lmps_df,
    period::Vector{Int64},
    bus_name::AbstractArray=["bus5"],
    node::String="bus5",
)
    lmps_df = sort(lmps_df)
    aux_period=[]
    for t in period
        aux_period=vcat(aux_period, DateTime(initial_time)+Hour(t-1))
    end
    indices = []
    max_element=0
    min_element=0
    data = Array{Any}(nothing, (length(period), length(bus_name) + 1, length(lmps_df), length(lmps_df[collect(keys(lmps_df))[1]]))) 
    for (i, v) in enumerate(keys(lmps_df))
        for (l, k) in enumerate(keys(lmps_df[collect(keys(lmps_df))[1]]))
            for t in 1:length(period)
                data[t, 1, i, l] = aux_period[t]
                c = 2
                for j in bus_name
                    prices_hour = lmps_df[v][k][
                        aux_period[t] .<= lmps_df[v][k].DateTime .< 
                        aux_period[t]+ Hour(1),
                        j,     
                    ]

                    data[t, c, i, l] = sum(prices_hour; dims=1)[1]

                    if data[t,c,i,l]>max_element && data[t,c,i,l]<1e3
                        max_element=data[t,c,i,l]
                    elseif data[t,c,i,l]<min_element && data[t,c,i,l]>-1e3
                        min_element=data[t,c,i,l]
                    end

                    c = c + 1
                end
            end
        end
        indices = vcat(indices, v)
    end
    palette = :Dark2_8
    c = 1
    for (l, k) in enumerate(keys(lmps_df[collect(keys(lmps_df))[1]])) #different Problems
        for b in 1:length(bus_name)
            for t in 1:length(period)
                if length(lmps_df[collect(keys(lmps_df))[1]])>1 && k=="RT" 
                    if c == 1
                        plot(
                            indices,
                            data[t, b + 1, :, l];
                            label="hour:" * string(period[t] - 1) * "- " * string(bus_name[b])* "- " * k ,
                            legend=:outertopright,
                            linestyle=:dash,
                            palette=palette,
                        )
                    else
                        plot!(
                            indices,
                            data[t, b + 1, :, l];
                            label="hour:" * string(period[t] - 1) * "- " * string(bus_name[b])* "- " * k,
                            legend=:outertopright,
                            linestyle=:dash,
                            palette=palette,
                        )
                    end
                else
                    if c == 1
                        plot(
                            indices,
                            data[t, b + 1, :, l];
                            label="hour:" * string(period[t] - 1) * "- " * string(bus_name[b])* "- " * k ,
                            legend=:outertopright,
                            palette=palette,
                        )
                    else
                        plot!(
                            indices,
                            data[t, b + 1, :, l];
                            label="hour:" * string(period[t] - 1) * "- " * string(bus_name[b])* "- " * k,
                            legend=:outertopright,
                            palette=palette,
                        )
                    end

                end
                c = c + 1
            end
        end
    end

    return plot!(;
        title="Price per Virtual Bid on " * node,
        ylabel="Prices (\$/MWh)",
        xlabel="Bid offers (p.u.)",
        ylims = (min_element*0.9,max_element*1.1)
    )
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
function plot_revenue_curves(
    lmps_df,
    results_df::Dict{Any,Any},
    market_simulator,
    period::Vector{Int64},
    generator_name::String,
)
    lmps_df = sort(lmps_df)
    gen = get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_name = get_name(get_bus(gen))

    indices = []
    aux_period=[]
    min_element=0
    max_element=0
    for t in period
        aux_period=vcat(aux_period, DateTime(initial_time)+Hour(t-1))
    end
    data = Array{Any}(nothing, (length(period), 2, length(lmps_df)))
    price=zeros(length(keys(lmps_df[collect(keys(lmps_df))[1]])))
    price=Dict()
    for k in keys(lmps_df[collect(keys(lmps_df))[1]])
        price[k]=0
    end
    
    for (i, v) in enumerate(keys(lmps_df))
        for t in 1:length(period)
            data[t, 1, i] = aux_period[t]
            variable_results = read_realized_variables(
                    get_problem_results(results_df[v]["ED"], "UC"); names=[:P__ThermalStandard]
                )
            generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
            virtual_gen = generator_data[1][!, generator_name][[period[t]]][1] 

            for k in (keys(lmps_df[collect(keys(lmps_df))[1]])) #Problems    
                prices_hour = lmps_df[v][k][
                    aux_period[t] .<= lmps_df[v][k].DateTime .< 
                    aux_period[t]+ Hour(1),
                    node,     
                ]
                price[k] = sum(prices_hour; dims=1)[1]
            end
            data[t, 2, i] = (price["DA"]-price["RT"])*virtual_gen 

            if data[t,2,i]>max_element && data[t,2,i]<1e3
                max_element=data[t,2,i]
            elseif data[t,2,i]<min_element && data[t,2,i]>-1e3
                min_element=data[t,2,i]
            end

        end
        indices = vcat(indices, v)
    end
    palette= :Dark2_8
    c = 1
    for t in 1:length(period)
        if c == 1
            plot(
                indices,
                data[t, 2, :];
                label="hour:" * string(period[t] - 1),
                legend=:outertopright,
                palette=palette,
            )
        else
            plot!(
                indices,
                data[t, 2, :];
                label="hour:" * string(period[t] - 1),
                legend=:outertopright,
                palette=palette,
            )
        end
        c = c + 1
    end

    return plot!(;
        title="Virtual Revenue per Offer on " * bus_name,
        ylabel="Revenue (\$/MWh)",
        xlabel="Bid offers (p.u)",
        ylims = (min_element*0.9,max_element*1.1)
    )
    #TODO: Change x axis 

end

"""
    plot_revenue_curves_renewable(
        lmps_df::Dict{Any, Any}, results_df::Dict{Any, Any}, market_simulator::MarketSimulator, period::bids{Float}, 
        bus_name::AbstractArray=["bus5"],
    )

Function to plot the revenue curve for the the renewable generators. 
The 'generator_name' defines which is the generator that we want to plot it's results
and 'bids' controls which bids we want to include in the plot.
"""
function plot_revenue_curves_renewable(
    lmps_df,
    results_df::Dict{Any,Any},
    market_simulator,
    bids::Vector{Float64},
    generator_name::String,
)
    gen = get_component(RenewableDispatch, market_simulator.system_uc, generator_name)
    bus_name = get_name(get_bus(gen))
    min_element=0
    max_element=0
    data = Array{Any}(nothing, (24, length(bids)+1))
    data[:,1] = lmps_df[first(keys(lmps_df))]["DA"][!,"DateTime"]
    for (j,q) in enumerate(bids) 
        variable_results = read_realized_variables(
                    get_problem_results(results_df[q]["ED"], "UC"); names=[:P__RenewableDispatch]
                )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_da = generator_data[1][!, generator_name] 

        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["RT"], "RT"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_rt_aux = generator_data[1][!, generator_name] 
        gen_rt=[]
        
        for i=1:(round(Int,length(gen_rt_aux)/12)) #TODO: Change to horizon
            gen_rt = vcat(gen_rt,[sum(gen_rt_aux[(1+12*(i-1)):(12+12*(i-1))])/12])
        end
        price=Dict()
        for k in (keys(lmps_df[collect(keys(lmps_df))[1]])) #Problems 
            for t =1:24 
                prices_hour = lmps_df[q][k][
                    lmps_df[first(keys(lmps_df))]["DA"][!,"DateTime"][t] .<= lmps_df[q][k].DateTime .< 
                    lmps_df[first(keys(lmps_df))]["DA"][!,"DateTime"][t]+ Hour(1),
                    bus_name,     
                ]
                if t == 1
                    price[k] = sum(prices_hour)
                else
                    price[k] = vcat(price[k],sum(prices_hour))
                end
            end
        end
          
        data[:,j+1] = gen_da.*price["DA"] + (gen_rt - gen_da).*price["RT"]

        if maximum(data[:,j+1])>max_element && minimum(data[:,j+1])<1e3
            max_element=maximum(data[:,j+1])
        elseif minimum(data[:,j+1])<min_element && minimum(data[:,j+1])>-1e3
            min_element=minimum(data[:,j+1])
        end
    end
    palette = :Dark2_8

    c = 1
    for (j,q) in enumerate(bids)
        if c == 1
            plot(
                0:length(data[:,1])-1,
                data[:, j+1];
                label="bid: " * string(q),
                legend=:outertopright,
                palette=palette,
            )
        else
            plot!(
                0:length(data[:,1])-1,
                data[:, j+1];
                label="bid: " * string(q),
                legend=:outertopright,
                palette=palette,
            )
        end
        c = c + 1
    end

    return plot!(;
        title= generator_name *" Revenue per Virtual Offer on " * node,
        ylabel="Revenue (\$/MWh)",
        xlabel="Period",
        ylims = (min_element*0.9,max_element*1.1)
    )
    #TODO: Change x axis 

end

"""
    plot_revenue_curves_renewable(
        lmps_df::Dict{Any, Any}, results_df::Dict{Any, Any}, market_simulator::MarketSimulator, period::bids{Float}, 
        bus_name::AbstractArray=["bus5"],
    )

Function to plot the revenue curve for the the renewable generators. 
The 'generator_name' defines which is the generator that we want to plot it's results
and 'bids' controls which bids we want to include in the plot.
"""
function plot_revenue_curves_renewable_plus_virtual(
    lmps_df,
    results_df::Dict{Any,Any},
    market_simulator,
    bids::Vector{Float64},
    renewable_gen::String,
    virtual_gen::String,
)
    gen_r = get_component(RenewableDispatch, market_simulator.system_uc, renewable_gen)
    bus_r = get_name(get_bus(gen_r))
    gen_v = get_component(ThermalStandard, market_simulator.system_uc, virtual_gen)
    bus_v = get_name(get_bus(gen_v))
    min_element=0
    max_element=0
    data = Array{Any}(nothing, (24, length(bids)+1,3))
    data[:,1,1] = lmps_df[first(keys(lmps_df))]["DA"][!,"DateTime"]
    data[:,1,2] = data[:,1,1]
    data[:,1,3] = data[:,1,1]
    for (j,q) in enumerate(bids) 
        variable_results = read_realized_variables(
                    get_problem_results(results_df[q]["ED"], "UC"); names=[:P__RenewableDispatch]
                )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_da_r = generator_data[1][!, renewable_gen] 

        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["RT"], "RT"); names=[:P__RenewableDispatch]
        )
        generator_data = getindex.(Ref(variable_results), [:P__RenewableDispatch])
        gen_rt_aux = generator_data[1][!, renewable_gen] 
        gen_rt_r=[]
        
        for i=1:(round(Int,length(gen_rt_aux)/12)) #TODO: Change to horizon
            gen_rt_r = vcat(gen_rt_r,[sum(gen_rt_aux[(1+12*(i-1)):(12+12*(i-1))])/12])
        end

        variable_results = read_realized_variables(
            get_problem_results(results_df[q]["ED"], "UC"); names=[:P__ThermalStandard]
        )
        generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
        gen_da_v = generator_data[1][!, generator_name]

        price=Dict()
        for k in (keys(lmps_df[collect(keys(lmps_df))[1]])) #Problems 
            for t =1:24 
                prices_hour = lmps_df[q][k][
                    lmps_df[first(keys(lmps_df))]["DA"][!,"DateTime"][t] .<= lmps_df[q][k].DateTime .< 
                    lmps_df[first(keys(lmps_df))]["DA"][!,"DateTime"][t]+ Hour(1),
                    :,     
                ]
                if t == 1
                    prices_hour[!,"DateTime"].=lmps_df[first(keys(lmps_df))]["DA"][!,"DateTime"][t]
                    price[k] = combine(groupby(prices_hour, :DateTime), names(prices_hour, Not(:DateTime)) .=> sum, renamecols=false)
                else
                    prices_hour[!,"DateTime"].=lmps_df[first(keys(lmps_df))]["DA"][!,"DateTime"][t]
                    price[k] = vcat(price[k],combine(groupby(prices_hour, :DateTime), names(prices_hour, Not(:DateTime)) .=> sum, renamecols=false))
                end
            end
        end
          
        data[:,j+1,1] = gen_da_r.*price["DA"][!,bus_r] + (gen_rt_r - gen_da_r).*price["RT"][!,bus_r]
        data[:,j+1,2] = gen_da_v.*(price["DA"][!,bus_v] - price["RT"][!,bus_v])

        if maximum(data[:,j+1,1]+data[:,j+1,2])>max_element && minimum(data[:,j+1,1]+data[:,j+1,2])<1e3
            max_element=maximum(data[:,j+1,1]+data[:,j+1,2])
        elseif minimum(data[:,j+1,1]+data[:,j+1,2])<min_element && minimum(data[:,j+1,1]+data[:,j+1,2])>-1e3
            min_element=minimum(data[:,j+1,1]+data[:,j+1,2])
        end
    end
    data[:,2:length(bids)+1,3]=data[:,2:length(bids)+1,1]+data[:,2:length(bids)+1,2]
    palette = :Dark2_8

    c = 1
    for (j,q) in enumerate(bids)
        if c == 1
            plot(
                0:length(data[:,1,1])-1,
                data[:, j+1,1];
                label="bid: " * string(q),
                legend=:outertopright,
                palette=palette,
            )
        else
            plot!(
                0:length(data[:,1,1])-1,
                data[:, j+1,1];
                label="bid: " * string(q),
                legend=:outertopright,
                palette=palette,
            )
        end
        c = c + 1
    end

    return plot!(;
        title= generator_name *" Revenue per Virtual Offer on " * node,
        ylabel="Revenue (\$/MWh)",
        xlabel="Period",
        ylims = (min_element*0.9,max_element*1.1)
    )
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

function plot_generation_curves(
    lmps_df,
    results_df::Dict{Any,Any},
    market_simulator,
    period::Vector{Int64},
    generator_name::String,
)
    lmps_df = sort(lmps_df)
    gen = get_component(ThermalStandard, market_simulator.system_uc, generator_name)
    bus_name = get_name(get_bus(gen))
    aux_period=[]
    for t in period
        aux_period=vcat(aux_period, DateTime(initial_time)+Hour(t-1))
    end

    indices = []
    data = Array{Any}(nothing, (length(period), 2, length(lmps_df)))
    for (i, v) in enumerate(keys(lmps_df))
        for t in 1:length(period)
            data[t, 1, i] = aux_period[t]
            variable_results = read_realized_variables(
                get_problem_results(results_df[v]["ED"], "UC"); names=[:P__ThermalStandard]
            )
            generator_data = getindex.(Ref(variable_results), [:P__ThermalStandard])
            virtual_gen = generator_data[1][!, generator_name][[period[t]]][1]
            data[t, 2, i] = virtual_gen #arrumar 
        end
        indices = vcat(indices, v)
    end
    palette= :Dark2_8

    c = 1
    for t in 1:length(period)
        if c == 1
            plot(
                indices,
                data[t, 2, :];
                label="hour:" * string(period[t] - 1),
                legend=:outertopright,
                palette = palette,
            )
        else
            plot!(
                indices,
                data[t, 2, :];
                label="hour:" * string(period[t] - 1),
                legend=:outertopright,
                palette = palette,
            )
        end
        c = c + 1
    end

    return plot!(;
        title=generator_name * " generation per Offer on " * bus_name,
        ylabel="Generation(MWh)",
        xlabel="Bid offers (p.u)",
    )
end
