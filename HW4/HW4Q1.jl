#=
 = Including reusable code
=#
include("HW4.jl")

#=
 = Part a - Method of Characteristics
=#
# Nothing to be done here

#=
 = Part b - Solution and Plotting
=#
f(x,t,a) = exp(-(x*exp(-t) - 1)^2) + a*x*(1-exp(-t))
save_plot(plot([(x) -> f(x,0.1,0),
		(x) -> f(x,2,0),
		(x) -> f(x,0.1,0.05),
		(x) -> f(x,2,0.05)],
		 -20, 20, color = ["a=0, t=0.1", "a=0, t=2", "a=0.05, t=0.1", "a=0.05, t=2"],
		 Guide.colorkey("Parameters"),
		 Guide.title("Solutions"),
		 Guide.YLabel("f"),
		 Guide.XLabel("x")), "q1")