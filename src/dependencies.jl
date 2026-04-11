struct DependencyNode
    cname::Symbol
    level::Int
    parent::Union{Nothing, Symbol}
    childs::Vector{DependencyNode}
    DependencyNode(cname::Symbol, level::Int, parent) = new(cname, level, parent, Vector{DependencyNode}())
end


function deptree(model::Model)
    function deptree(model, cname::Symbol, level::Int, parent::Union{Nothing, Symbol})
        out = DependencyNode(cname, level, parent)
        for d in dependencies(model, cname)
            push!(out.childs, deptree(model, d, level+1, cname))
        end
        return out
    end

    # Identify parent for all comps
    parent = OrderedDict{Symbol, Symbol}()
    for (cname, comp) in model.comps
        for d in GModelFit.dependencies(model, cname)
            @assert !haskey(parent, d) "Component $d has two parent nodes: $(parent[d]) and $cname"
            parent[d] = cname
        end
    end

    # Ensure no circular dependency is present by checking all parent
    # nodes of a given component to be different from the component
    # itself.  Also collect components with no parent.
    comps_with_no_parent = Vector{Symbol}()
    for cname in keys(model.comps)
        if haskey(parent, cname)
            p = parent[cname]
            @assert p != cname "Component $cname depends on itself"
            while haskey(parent, p)
                p = parent[p]
                if cname == p
                    display(parent)
                    error("Circular dependency detected for component $cname")
                end
            end
        else
            push!(comps_with_no_parent, cname)
        end
    end

    # Neglect components with no parent and no dependencies
    while length(comps_with_no_parent) > 1
        if length(dependencies(model, comps_with_no_parent[1])) == 0
            deleteat!(comps_with_no_parent, 1)
        end
    end

    # Main component
    if isnothing(model.maincomp)
        maincomp = comps_with_no_parent[end]
    else
        maincomp = model.maincomp
    end

    return deptree(model, maincomp, 1, nothing)
end


function flatten(node::DependencyNode)
    out = [node]
    for child in node.childs
        append!(out, flatten(child))
    end
    return out
end


function dependencies(model::Model, cname::Symbol)
    output = Vector{Symbol}()
    comp = model.comps[cname]
    for d in dependencies(comp)
        @assert haskey(model.comps, d)  ||  isa(comp, FComp)  "No component has name $d"
        haskey(model.comps, d)  &&  push!(output, d)
    end
    return output
end
