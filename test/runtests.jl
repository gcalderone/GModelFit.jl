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


# ====================================================================
function simulate_measures(domain::AbstractDomain{N}, val::AbstractArray{T, N}, noise=0.05) where {T <: Real, N}
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


# Tie two parameters
model[:l2].norm.patch = :l1
res = fit!(model, data)

# Patch one parameter to another via a λ function
model[:l2].norm.fixed = false
model[:l2].norm.patch = @λ (v, m) -> v + m[:l1].norm
res = fit!(model, data)



# ====================================================================
x = 0:0.05:6
model1 = Model(Domain(x),
               :l1  => GFit.Gaussian(1, 2, 0.2),
               :l2  => GFit.Gaussian(1, 3, 0.5),
               :bkg => GFit.OffsetSlope(0.5, 1, 0.1),
               :main => SumReducer(:l1, :l2, :bkg));

model2 = Model(Domain(x),
               :l1  => GFit.Gaussian(0.8, 2.1, 0.1),
               :l2  => GFit.Gaussian(1.2, 2.5, 0.4),
               :bkg => GFit.OffsetSlope(0.5, 1, 0.1),
               :main => SumReducer(:l1, :l2, :bkg));

y1 = model1();
y2 = model2();
data1 = simulate_measures(Domain(x), y1)
data2 = simulate_measures(Domain(x), y2)

model = MultiModel(model1, model2)
res = fit!(model, [data1, data2])




# Global model
model = Model([pred1, pred2]);
freeze(model, :bkg);

@gp "set grid" tit="Multi-epoch data fitting" linetypes(:tab10, lw=1.5) "set bars 0"
@gp :- x model[1]() "w l t 'True model (T1)'"
@gp :- x model[2]() "w l t 'True model (T2)'"

# Simulate data measurements
data1 = Measures(model[1]() .+ noise .* randn(length(x)), noise);
data2 = Measures(model[2]() .+ noise .* randn(length(x)), noise);


@gp :- x data1.val data1.unc "w yerr t 'Data (T1)'"
@gp :- x data2.val data2.unc "w yerr t 'Data (T2)'"

model[1][:l2].center.fixed = true
patch!(model) do m
    m[1][:l2].center = m[2][:l2].center
end

@time res = fit!(model, [data1, data2])


res.cost/res.dof

@gp :- x model[1]() "w l t 'Best fit model (T1)'"
@gp :- x model[2]() "w l t 'Best fit model (T2)'"

thaw(model, :bkg);
@time res = fit!(model, [data1, data2]);


viewer(model, [data1, data2], res, selcomps=[:bkg, :l1, :l2])
