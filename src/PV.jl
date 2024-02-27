
module PV

using DataStructures

import Base.keys, Base.getindex, Base.setindex!,
       Base.propertynames, Base.getproperty, Base.setproperty!,
       Base.push!, Base.empty!, Base.iterate

export PVComp, PVModel, items

struct PVComp{T}
    pnames::Vector{Symbol}
    indices::Vector{Int}
    data::Vector{T}
end

struct PVModel{T}
    comps::OrderedDict{Symbol, PVComp{T}}
    indices::Vector{Int}
    data::Vector{T}
end


PVComp(parent::PVModel{T}) where T = PVComp{T}(Vector{Symbol}(), Vector{Int}(), parent.data)
propertynames(comp::PVComp) = getfield(comp, :pnames)

function index(comp::PVComp, pname::Symbol)
    i = findfirst(getfield(comp, :pnames) .== pname)
    @assert !isnothing(i) "Unknwon parameter name: $pname"
    return getfield(comp, :indices)[i]
end
getindex(comp::PVComp, pname::Symbol) = getfield(comp, :data)[index(comp, pname)]
setindex!(comp::PVComp, value, pname::Symbol) = setindex!(getfield(comp, :data),
                                                          value, index(comp, pname))
getproperty( comp::PVComp, pname::Symbol) = getindex(comp, pname)
setproperty!(comp::PVComp, pname::Symbol, value) = setindex!(comp, value, pname)

items(comp::PVComp) = view(getfield(comp, :data), getfield(comp, :indices))

function iterate(comp::PVComp, i=1)
    (i > length(getfield(comp, :pnames)))  &&  return nothing
    return (getfield(comp, :pnames)[i] => getfield(comp, :data)[getfield(comp, :indices)[i]], i+=1)
end



PVModel{T}() where T = PVModel{T}(OrderedDict{Symbol, PVComp{T}}(), Vector{Int}(), Vector{T}())

function empty!(pv::PVModel)
    empty!(pv.comps)
    empty!(pv.indices)
    empty!(pv.data)
end

keys(pv::PVModel) = collect(keys(pv.comps))

function push!(pv::PVModel{T}, cname::Symbol, pname::Symbol, value::T) where T
    comp = pv[cname]
    if pname in getfield(comp, :pnames)
        comp[pname] = value
    else
        i = length(getfield(comp, :data)) + 1
        push!(getfield(comp, :pnames), pname)
        push!(getfield(comp, :indices), i)
        push!(getfield(comp, :data), value)
    end
    # Collect indices from components
    empty!(pv.indices)
    for (cname, comp) in pv.comps
        append!(pv.indices, getfield(comp, :indices))
    end
    return value
end

function getindex(pv::PVModel{T}, cname::Symbol) where T
    (cname in keys(pv.comps))  ||  (pv.comps[cname] = PVComp(pv))
    return pv.comps[cname]
end

items(pv::PVModel) = view(pv.data, pv.indices)
iterate(pv::PVModel, state...) = iterate(pv.comps, state...)

end
