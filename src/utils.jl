using Dierckx, QuadGK


function int_tabulated(λ, f; int_k=3)
    @assert issorted(λ)
    @assert all(isfinite.(λ))
    @assert all(isfinite.(f))
    itp = Spline1D(λ, f, k=int_k, bc="error")
    return quadgk(itp, λ[1], λ[end])
end
