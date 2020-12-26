struct CDomain <: AbstractComponent
    dim::Int
    CDomain(dim::Int=1) = new(dim)
end

compeval_array(comp::CDomain, domain::AbstractDomain) = deepcopy(domain[comp.dim])

evaluate(buffer, comp::CDomain, domain::AbstractDomain) = nothing
