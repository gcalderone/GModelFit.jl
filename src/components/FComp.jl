
struct FComp <: GFit.AbstractComponent
    params::Vector{Parameter}
    f::Function
    FComp(npar::Int, funct::Function) = new([Parameter(1.) for i in 1:npar], funct)
end

function evaluate!(buffer::Vector{Float64}, comp::FComp, x::AbstractDomain,
                   params::Vararg{Float64})
    buffer .= getfield(comp, :f)([params...])
end
