abstract type AbstractFitData end

struct FitData <: AbstractFitData
    timestamp::DateTime
    model::Model
    data::Vector{Float64}
    unc::Vector{Float64}
    resid::Vector{Float64}
    ifree::Vector{Int}
    dof::Int
    progress::ProgressUnknown

    function FitData(model::Model, data::Measures{N}) where N
        evaluate(model)
        data1d = flatten(data)
        resid = (model() .- data1d.val) ./ data1d.unc
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
        @assert length(ifree) > 0 "No free parameter in the model"

        dof = length(resid) - length(ifree)
        prog = ProgressUnknown("Model (dof=$dof) evaluations:", dt=0.5, showspeed=true)
        return new(now(), model, data1d.val, data1d.unc, resid, ifree, dof, prog)
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



# TODO: MultiFitData
