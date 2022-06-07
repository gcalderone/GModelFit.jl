using Random, Test, GFit

hv = GFit.HashVector{Int}()
hv.a = 1
hv.b = 2
hv.b = 3

@test length(hv) == 2
@test all(propertynames(hv) .== [:a, :b])

@test hv.a == 1
@test hv.b == 3

hv.a = 10
@test hv.a == 10

for (key, val) in hv
    println(key, "  ", val)
end

@test hv[:a] == 10
@test hv[:b] == 3

hv[:a] = 88
hv[:b] = 99

@test hv[:a] == 88
@test hv[:b] == 99


hhv = GFit.HashHashVector{Int}()
hhv[:comp1][:a] = 10
hhv[:comp1][:b] = 20
hhv[:comp2][:a] = 30


@test hhv[:comp1].a == 10
@test hhv[:comp1].b == 20
@test hhv[:comp2].a == 30

@test hhv[:comp1][:a] == 10
@test hhv[:comp1][:b] == 20
@test hhv[:comp2][:a] == 30


for (key1, hv) in hhv
    for (key, val) in hv
        println(key1, "  ", key, "  ", val)
    end
end


function simulate_measures(domain::AbstractDomain{N}, val::AbstractArray{T, N}, noise=0.1) where {T <: Real, N}
    ee = extrema(val)
    range = ee[2] - ee[1]
    Measures(domain, val .+ noise .* range .* randn(length(val)), noise .* range)
end


# ====================================================================
x = 0:0.1:5
y = x.^2 .+ x .+ 5;
data = Measures(Domain(x), y .+ randn(length(x)), 1.)
model = Model(Domain(x), :a2 => 1, :a1 => 1, :a0 => 5)
model[:parabola] = @λ (x, a2, a1, a0) -> @. (a2 * x^2  +  a1 * x  + a0)
res = fit!(model, data)
# @gp x y "w l t 'True model'" x data.val data.unc "w yerr t 'Data'" x model() "w l t 'Best fit'"


# ====================================================================
x = 0:0.1:5
y = x.^2 .+ x .+ 5;
data = simulate_measures(Domain(x), y)
model = Model(Domain(x), @λ (x, a2=1, a1=1, a0=5) -> @. (a2 * x^2  +  a1 * x  + a0))
res = fit!(model, data)


# ====================================================================
f = @λ (x, p1=1, p2=1.e-3, p3=1e-6, p4=4, p5=5) ->
@. (p1 + p2 * x + p3 * x^2 + p4 * sin(p5 * x))  *  cos(x)

x = 1.:50:10000
y = f(x);
data = simulate_measures(Domain(x), y)
model = Model(Domain(x), f)
res = fit!(model, data)


# ====================================================================
f1 = @λ (x, p1=1, p2=1e-3, p3=1e-6) -> @.  p1  +  p2 * x  +  p3 * x^2
f2 = @λ (x, p4=4, p5=5) -> @. p4 * sin(p5 * x)
f3 = @λ (x) -> cos.(x)

x = 1.:50:10000
y = (f1.(x) .+ f2.(x)) .* f3.(x);
data = simulate_measures(Domain(x), y)

model = Model(Domain(x),
              :f1 => f1, 
              :f2 => f2,
              :f3 => f3,
              :main => @λ (x, f1, f2, f3) -> (f1 .+ f2) .* f3)
res = fit!(model, data)

# Same results with
model = Model(Domain(x), :f1 => f1)
model[:f2] = f2
model[:f3] = f3
model[:main] = @λ (x, f1, f2, f3) -> (f1 .+ f2) .* f3
res = fit!(model, data)



# ====================================================================
x = 0:0.05:6
model = Model(Domain(x),
              :l1  => GFit.Gaussian(1, 2, 0.2),
              :l2  => GFit.Gaussian(1, 3, 0.5),
              :bkg => GFit.OffsetSlope(0.5, 1, 0.1),
              :main => SumReducer(:l1, :l2, :bkg));
y = model();
data = simulate_measures(Domain(x), y)
res = fit!(model, data)
