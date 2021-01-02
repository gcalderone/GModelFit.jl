struct CDomain <: AbstractComponent
    dim::Int
    CDomain(dim::Int=1) = new(dim)
end

prepare!(comp::CDomain, domain::AbstractDomain) = domain[comp.dim]

evaluate!(buffer, comp::CDomain, domain::AbstractDomain) = nothing
