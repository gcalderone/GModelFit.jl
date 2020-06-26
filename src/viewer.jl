function dump(model::Model, args...; kw...)
    io = IOBuffer()
    dump(io, model, args...; kw...)
    return String(take!(io))
end

function dump(filename::String, model::Model, args...; kw...)
    io = open(filename, "w")
    dump(io, model, args...; kw...)
    close(io)
    return filename
end

function viewermeta(model::Model, data::AbstractData;
                    names=Vector{Symbol}(), index=1)
    default_meta = OrderedDict()
    default_meta["Label"] = "Label"
    default_meta["visible"] = true
    default_meta["color"] = "auto"

    meta = OrderedDict()
    meta["Title"] = "My Title"
    meta["Label_x"] = ""
    meta["Label_y"] = ""

    meta["data"] = deepcopy(default_meta)
    meta["data"]["Label"] = "Empirical data"
    meta["data"]["color"] = "black"

    pred = model.preds[index]
    if length(names) == 0
        names = collect(keys(pred.reducers))
    end
    meta["model"] = OrderedDict()
    for i in 1:length(names)
        n = string(names[i])
        meta["model"][n] = deepcopy(default_meta)
        meta["model"][n]["Label"] = "Model: $n"
    end

    meta["residuals"] = deepcopy(default_meta)
    meta["residuals"]["Label"] = "Normalized residuals"
    meta["residuals"]["color"] = "gray"

    return meta        
end

function dump(io::IO, model::Model, data::AbstractData;
              names=Vector{Symbol}(), index=1,
              format=:JSON,
              meta::AbstractDict=Dict())
    pred = model.preds[index]

    root = OrderedDict()
    root["meta"] = meta

    root["data"] = OrderedDict()
    root["data"]["x"] = domain(pred)
    root["data"]["y"] = data.val
    root["data"]["err"] = data.unc

    if length(names) == 0
        names = collect(keys(pred.reducers))
    end
    root["model"] = OrderedDict()
    for i in 1:length(names)
        n = string(names[i])
        root["model"][n] = OrderedDict()
        root["model"][n]["x"] = domain(pred)
        root["model"][n]["y"] = pred.revals[names[i]]
    end

    root["residuals"] = OrderedDict()
    root["residuals"]["x"] = domain(pred)
    root["residuals"]["y"] = (data.val .- pred()) ./ data.unc

    JSON.print(io, root)
end
