struct SumReducer <: AbstractComponent
    list::Vector{Symbol}
end

deps(comp::SumReducer) = comp.list

function evaluate!(buffer::Vector{Float64}, comp::SumReducer, domain::AbstractDomain,
                   args)
    buffer .= 0.
    for i in 1:length(args)
        buffer .+= args[i]
    end
end
