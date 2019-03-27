using LsqFit            
@. model(x, p) = p[1]*exp(-x*p[2])            
xdata = range(0, stop=10, length=1000)            
ydata = model(xdata, [1.0 2.0]) + 0.01*randn(length(xdata))            
p0 = [0.5, 0.5]            
fit = curve_fit(model, xdata, ydata, p0)            
fit.param # gives the fitted parameters p[1] and p[2]            
plot(xdata, ydata)
y1(x) = fit.param[1]*exp.(-x*fit.param[2])
plot!(xdata, y1(xdata))