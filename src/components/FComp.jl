struct FComp <: GFit.AbstractComponent
    funct::Function
    deps::Vector{Symbol}
    params::Vector{Parameter}

    FComp(funct::Function, guess::Vector{T}) where T <: Number =
        new(funct, Symbol[], Parameter.(guess))
    FComp(funct::Function, deps::Vector{Symbol}, guess::Vector{T}) where T <: Number =
        new(funct, deps, Parameter.(guess))
end

dependencies(comp::λComp) = getfield(comp, :deps)


# We need to implement two evaluate! methods, with/without deps argument respectively
function evaluate!(buffer::Vector{Float64}, comp::FComp, x::AbstractDomain,
                   params::Vararg{Float64})
    buffer .= getfield(comp, :funct)([params...])
end

function evaluate!(buffer::Vector{Float64}, comp::FComp, domain::AbstractDomain,
                   deps, params...)
    buffer .= getfield(comp, :func)(deps..., [params...])
end
