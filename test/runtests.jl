using Random, Test, GModelFit, GModelFit.PV

# Test Model
mm = PVModel{Float64}()
@assert typeof(mm) == PVModel{Float64}

mm[:comp]  # empty component
push!(mm, :comp3, :par, 3.1)
push!(mm, :comp1, :par, 1.1)
push!(mm, :comp1, :alt, 1.2)
mm[:comp3][:par] = 99

@assert getfield(mm, :data) == [99, 1.1, 1.2]

@assert mm[:comp1].par == 1.1
@assert mm[:comp1].alt == 1.2
@assert mm[:comp3].par == 99

@assert keys(mm) == [:comp, :comp3, :comp1]

@assert propertynames(mm[:comp])  == Symbol[]
@assert propertynames(mm[:comp1]) == [:par, :alt]
@assert propertynames(mm[:comp3]) == [:par]

@assert items(mm) == [99, 1.1, 1.2]
@assert items(mm[:comp]) == Float64[]
@assert items(mm[:comp1]) == [1.1, 1.2]
@assert items(mm[:comp3]) == [99]

for (cname, comp) in mm
    println(cname)
    for (pname, par) in comp
        println("  ", pname, ": ", par)
        @assert getproperty(mm[cname], pname) == par
    end
end



# ====================================================================
x    = [0.1, 1.1, 2.1, 3.1, 4.1]
meas = [6.29, 7.27, 10.41, 18.67, 25.3]
unc  = [1.1, 1.1, 1.1, 1.2, 1.2]

domain = Domain(x)
data = Measures(domain, meas, unc)
model = Model(@fd (x, a2=1, a1=1, a0=5) -> (a2 .* x.^2  .+  a1 .* x  .+  a0))
res = fit(model, data)
# @gp data model


# ====================================================================
x = 0:0.1:5
model = Model(:parabola => @fd (x, a2=1, a1=1, a0=5) -> @. (a2 * x^2  +  a1 * x  + a0))
res = fit(model, GModelFit.mock(Measures, model, Domain(x), seed=1))
# @gp x y "w l t 'True model'" x values(data) uncerts(data) "w yerr t 'Data'" x model() "w l t 'Best fit'"


# ====================================================================
x = 0:0.1:5
model = Model(@fd (x, a2=1, a1=1, a0=5) -> @. (a2 * x^2  +  a1 * x  + a0))
res = fit(model, GModelFit.mock(Measures, model, Domain(x), seed=1))


# ====================================================================
f = @fd (x, p1=1, p2=1.e-3, p3=1e-6, p4=4, p5=5) ->
@. (p1 + p2 * x + p3 * x^2 + p4 * sin(p5 * x))  *  cos(x)
x = 1.:50:10000
model = Model(f)
res = fit(model, GModelFit.mock(Measures, model, Domain(x), seed=1))


# ====================================================================
f1 = @fd (x, p1=1, p2=1e-3, p3=1e-6) -> @.  p1  +  p2 * x  +  p3 * x^2
f2 = @fd (x, p4=4, p5=5) -> @. p4 * sin(p5 * x)
f3 = @fd (x) -> cos.(x)

x = 1.:50:10000
model = Model(:f1 => f1,
              :f2 => f2,
              :f3 => f3,
              :main => @fd (x, f1, f2, f3) -> (f1 .+ f2) .* f3)
res = fit(model, GModelFit.mock(Measures, model, Domain(x), seed=1))

# Same results with
model = Model(:f1 => f1)
model[:f2] = f2
model[:f3] = f3
model[:main] = @fd (x, f1, f2, f3) -> (f1 .+ f2) .* f3
res = fit(model, GModelFit.mock(Measures, model, Domain(x), seed=1))


# ====================================================================
x = 0:0.05:6
model = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
              :l2  => GModelFit.Gaussian(1, 3, 0.5),
              :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
              :main => SumReducer(:l1, :l2, :bkg));
res = fit(model, GModelFit.mock(Measures, model, Domain(x), seed=1))


# Tie two parameters
model[:l2].norm.patch = :l1
res = fit(model, GModelFit.mock(Measures, model, Domain(x), seed=1))

# Patch one parameter to another via a λ function
model[:l2].norm.patch = @fd (m, v) -> v + m[:l1].norm
res = fit(model, GModelFit.mock(Measures, model, Domain(x), seed=1))



# ====================================================================
x = 0:0.05:6
model1 = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
               :l2  => GModelFit.Gaussian(1, 3, 0.5),
               :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
               :main => SumReducer(:l1, :l2, :bkg));

model2 = Model(:l1  => GModelFit.Gaussian(0.8, 2.1, 0.1),
               :l2  => GModelFit.Gaussian(1.2, 2.5, 0.4),
               :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
               :main => SumReducer(:l1, :l2, :bkg));

models = [model1, model2]
freeze!(models[1], :bkg);
freeze!(models[2], :bkg);
data = GModelFit.mock(Measures, models, [Domain(x), Domain(x)], seed=1)
res = fit(models, data, minimizer=GModelFit.cmpfit())
# GModelFit.print_param_covariance(res, sort=true, select=["[2][l1].norm"])



thaw!(models[1], :bkg);
thaw!(models[2], :bkg);

models[2][:bkg].offset.mpatch = @fd m -> m[1][:bkg].offset
models[2][:bkg].slope.mpatch  = @fd m -> m[1][:bkg].slope


models[1][:l2].center.mpatch = @fd m -> m[2][:l2].center

@time res = fit(models, data)

#=
@gp x y1 "w l t 'True model'" x values(data1) uncerts(data1) "w yerr t 'Data'" x model1() "w l t 'Best fit'"
@gp x y2 "w l t 'True model'" x values(data2) uncerts(data2) "w yerr t 'Data'" x model2() "w l t 'Best fit'"
=#
