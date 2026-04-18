struct DependencyNode
    cname::Symbol
    level::Int
    parent::Union{Nothing, Symbol}
    childs::Vector{DependencyNode}
    DependencyNode(cname::Symbol, level::Int, parent) = new(cname, level, parent, DependencyNode[])
end


function deptree(model::Model)
    function build_tree(cname::Symbol, level::Int, p::Union{Nothing, Symbol})
        node = DependencyNode(cname, level, p)
        for d in dependencies(model, cname)
            push!(node.childs, build_tree(d, level + 1, cname))
        end
        return node
    end

    # Identify parent for all comps
    parent = OrderedDict{Symbol, Symbol}()
    for cname in keys(model.comps)
        for d in dependencies(model, cname)
            @assert !haskey(parent, d) "Component $d has two parent nodes: $(parent[d]) and $cname"
            parent[d] = cname
        end
    end

    # Ensure no circular dependency is present and collect root components
    comps_with_no_parent = Symbol[]
    for cname in keys(model.comps)
        if haskey(parent, cname)
            p = parent[cname]
            visited = Set{Symbol}([cname])
            while haskey(parent, p)
                @assert p ∉ visited "Circular dependency detected for component $p"
                push!(visited, p)
                p = parent[p]
            end
        else
            push!(comps_with_no_parent, cname)
        end
    end

    # Determine the main component
    if !isnothing(model.maincomp)
        maincomp = model.maincomp
    else
        # Consider as possible main comp the ones with no parent but having dependencies
        valid_roots = filter(c -> !isempty(dependencies(model, c)), comps_with_no_parent)
        if isempty(valid_roots)
            maincomp = comps_with_no_parent[end]
        else
            maincomp = valid_roots[end]
        end
    end

    return build_tree(maincomp, 1, nothing)
end


function flatten(node::DependencyNode)
    out = [node]
    for child in node.childs
        append!(out, flatten(child))
    end
    return out
end


function dependencies(model::Model, cname::Symbol)
    output = Symbol[]
    comp = model.comps[cname]
    for d in dependencies(comp)
        @assert haskey(model.comps, d)  ||  isa(comp, FComp)  "No component has name $d"
        haskey(model.comps, d)  &&  push!(output, d)
    end
    return output
end
