struct CDomain <: AbstractComponent
    dim::Int
    CDomain(dim::Int=1) = new(dim)
end

compeval_cdata(comp::CDomain, domain::AbstractDomain) = nothing
compeval_array(comp::CDomain, domain::AbstractDomain) = deepcopy(flatten(domain)[comp.dim])

evaluate(c::CompEval{CDomain, T}) where T <: AbstractLinearDomain = nothing
