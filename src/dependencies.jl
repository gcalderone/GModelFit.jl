find_sequence(lh::LikelihoodSnapshot{GT}, gkey::GT) where {GT} = getgroups(lh)[gkey].seq
function find_sequence(model::AbstractModel{GT}, gkey::GT) where GT
    # Identify "active" components and their dependencies
    active = OrderedDict{Symbol, Vector{Symbol}}()
    for ckey in keys(getcomps(model))
        (group_key(model, ckey) == gkey)  ||  continue
        cname = comp_name(model, ckey)
        for dep in dependencies(getcomps(model)[ckey])
            haskey(active, cname)             ||  (active[cname] = Vector{Symbol}())
            if haskey(getcomps(model), (gkey..., dep))
                haskey(active, dep)           ||  (active[dep]   = Vector{Symbol}())
                push!(active[cname], dep)
            else
                error("Missing dependency: component $ckey depends on :$dep, but the latter is not defined.")
            end
        end
    end
    if length(active) == 0  # there are no depdenencies, just take the last component
        return [comp_name(model, collect(keys(getcomps(model)))[end])]
    end

    sequence = Symbol[]
    while !isempty(active)
        lbefore = length(sequence)
        for (cname, deps) in active  # Search for components with no dependencies
            isempty(deps)  ||  continue
            push!(sequence, cname)   # push into output sequence
            delete!(active, cname)   # delete from dictionary
            for (k, v) in active     # ...also delete from other component dependencies
                filter!(d -> d != cname, v)
            end
            break
        end
        (lbefore == length(sequence))  &&  error("""Circular dependency detected among components: $(join(keys(active), ", ")) in group $(gkey)""")
    end
    return sequence
end
