using Gnuplot, Dates
Gnuplot.quitall()
mkpath("assets")
Gnuplot.options.term = "unknown"
empty!(Gnuplot.options.init)
push!( Gnuplot.options.init, linetypes(:Set1_5, lw=2.5, ps=1.5))
function saveas(file)
    Gnuplot.save(term="pngcairo size 550,350 fontscale 0.8", "assets/$(file).png")
    nothing
end


function comparedata(A, B)
    function iscomparable(a)
        tt = typeof(a)
        return ((tt <: Number)          ||
                (tt <: AbstractString)  ||
                (tt <: Enum)            ||
                (tt <: Dates.AbstractTime))
    end
    TA = typeof(A)
    TB = typeof(B)
    if TA != TB
        @warn "$TA != $TB"
    else
        if (TA <: AbstractVector)  ||
            (TA <: Tuple)
            @assert length(A) == length(B)
            for i in 1:length(A)
                a = getindex(A, i)
                b = getindex(B, i)
                if iscomparable(a)
                    if a != b
                        @warn TA i a b
                    end
                else
                    comparedata(a, b)
                end
            end
        elseif isstructtype(TA)
            for i in 1:nfields(A)
                a = getfield(A, i)
                b = getfield(B, i)
                if iscomparable(a)
                    if a !== b
                        if isa(a, Number)  &&  !isnan(a)  &&  !isnan(b)
                            @warn TA fieldname(TA, i) a b
                        end
                    end
                else
                    comparedata(a, b)
                end
            end
        else
            @warn "Unhandled type: $TA"
        end
    end
end


function dumpjson(file, args...)
    restored = GModelFit.deserialize(GModelFit.serialize("assets/$(file).json", args...))
    rm("assets/$(file).json")
    # comparedata(restored, [args...])
end
