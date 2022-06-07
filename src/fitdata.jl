abstract type AbstractFitData end

# ====================================================================
struct FitData <: AbstractFitData
    timestamp::DateTime
    progress::ProgressUnknown
    model::Model
    data::Vector{Float64}
    unc::Vector{Float64}
    resid::Vector{Float64}
    ifree::Vector{Int}
    dof::Int

    function FitData(model::Model, data::Measures{N}) where N
        evaluate(model)
        data1d = flatten(data)
        resid = fill(NaN, length.(data1d))
        ifree = Vector{Int}()
        i = 1
        for (cname, hv) in model.params
            for (pname, par) in hv
                if !par.fixed  &&  (model.cevals[cname].cfixed == 0)
                    push!(ifree, i)
                end
                i += 1
            end
        end
        nfree = length(ifree)
        @assert nfree > 0 "No free parameter in the model"
        dof = length(resid) - nfree
        prog = ProgressUnknown("Model (dof=$dof) evaluations:", dt=0.5, showspeed=true)
        return new(now(), prog, model, data1d.val, data1d.unc, resid, ifree, dof)
    end
end

free_params(fd::FitData) = internal_data(fd.model.params)[fd.ifree]
residuals(fd::FitData) = fd.resid

function evaluate!(fd::FitData, pvalues::Vector{Float64})
    internal_data(fd.model.pvalues)[fd.ifree] .= pvalues
    eval_step2(fd.model)
    eval_step3(fd.model)
    fd.resid .= (fd.model() .- fd.data) ./ fd.unc
    ProgressMeter.next!(fd.progress; showvalues=() -> begin
                        [(:fit_stat, sum(abs2, fd.resid) / fd.dof)]
                        end)
    return fd.resid
end

function finalize!(fd::FitData, best::Vector{Float64}, unc::Vector{Float64})
    @assert length(fd.ifree) == length(best) == length(unc)
    evaluate!(fd, best)
    all_unc = fill(NaN, length(internal_data(fd.model.params)))
    all_unc[fd.ifree] = unc
    eval_step4(fd.model, all_unc)
    ProgressMeter.finish!(fd.progress)
end

function error!(fd::FitData)
    eval_step4(fd.model)
    ProgressMeter.finish!(fd.progress)
end



# ====================================================================
struct MultiFitData <: AbstractFitData
    timestamp::DateTime
    progress::ProgressUnknown
    multi::MultiModel
    fds::Vector{FitData}
    resid::Vector{Float64}
    dof::Int

    function MultiFitData(multi::MultiModel, datasets::Vector{Measures{N}}) where N
        @assert length(multi) == length(datasets)
        evaluate(multi)
        fds = [FitData(multi[id], datasets[id]) for id in 1:length(multi)]
        resid = fill(NaN, sum(length.(getfield.(fds, :data))))
        nfree = sum(length.(getfield.(fds, :ifree)))
        @assert nfree > 0 "No free parameter in the model"
        dof = length(resid) - nfree
        prog = ProgressUnknown("Model (dof=$dof) evaluations:", dt=0.5, showspeed=true)
        return new(now(), prog, multi, fds, resid, dof)
    end
end

function free_params(fd::MultiFitData)
    out = Vector{Parameter}()
    for id in 1:length(fd.multi)
        append!(out, free_params(fd.fds[id]))
    end
    return out
end
residuals(fd::MultiFitData) = fd.resid

function evaluate!(fd::MultiFitData, pvalues::Vector{Float64})
    # We need to copy all parameter values before evaluate!, to ensure
    # all patch functions use the current parameter values
    i1 = 1
    for id in 1:length(fd.multi)
        nn = length(fd.fds[id].ifree)
        if nn > 0
            i2 = i1 + nn - 1
            internal_data(fd.fds[id].model.pvalues)[fd.fds[id].ifree] .= pvalues[i1:i2]
            i1 += nn
        end
    end
    # Now proceed to evaluation
    i1 = 1
    for id in 1:length(fd.multi)
        nn = length(fd.fds[id].ifree)
        if nn > 0
            i2 = i1 + nn - 1
            evaluate!(fd.fds[id], pvalues[i1:i2])
            i1 += nn
        end
    end
    # Populate resid vector
    i1 = 1
    for id in 1:length(fd.multi)
        nn = length(fd.fds[id].resid)
        if nn > 0
            i2 = i1 + nn - 1
            fd.resid[i1:i2] .= fd.fds[id].resid
            i1 += nn
        end
    end
    return fd.resid
end


function finalize!(fd::MultiFitData, best::Vector{Float64}, unc::Vector{Float64})
    i1 = 1
    for id in 1:length(fd.multi)
        nn = length(fd.fds[id].ifree)
        if nn > 0
            i2 = i1 + nn - 1
            finalize!(fd.fds[id], best[i1:i2], unc[i1:i2])
            i1 += nn
        end
    end
    ProgressMeter.finish!(fd.progress)
end


function error!(fd::MultiFitData)
    for id in 1:length(fd.multi)
        finalize!(fd.fds[id])
    end
    ProgressMeter.finish!(fd.progress)
end
