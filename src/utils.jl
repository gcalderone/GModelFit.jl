
function int_tabulated(x, y)
    @assert issorted(x)
    @assert all(isfinite.(x))
    @assert all(isfinite.(y))
    b = x[2:end] .- x[1:end-1]
    h = y[2:end] .+ y[1:end-1]
    return sum(b .* h) / 2
end
