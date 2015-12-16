#=
 = Including reusable code
=#
include("HW5.jl")

#=
 = Setup
=#
# Constant parameters of the problem
const N = 200
const epsilon = 1.0
const beta = 2 #beta = 1/kT = 1/(1/2)
const max_iter = 20000 #After some fiddling this seems to be a good range

#=
 = Part a - Setup
=#
particles = ones(N)
energy(particles) = sum(particles)

#=
 = Part b - Metropolis algorithm
=#

function metropolis(particles)
	part = floor(Int, rand()*N) + 1
	step = epsilon*(2*rand() - 1)
	if step + particles[part] < 0
		return
	elseif step <= 0
		particles[part] += step
	else
		if rand() < exp(-beta*step)
			particles[part] += step
		end
	end
end #metropolis

energies = zeros(max_iter)

#=
 = Part c - equil time
=#
const equil_time = 1000

for iter = 1:max_iter
	metropolis(particles)
	energies[iter] = energy(particles)
end

u_infinity = N/2
u_initial = N

model(t) = (u_infinity+(u_initial - u_infinity)*exp(-t/equil_time))

save_plot(plot(layer(x=1:max_iter, y=energies, Geom.line),
		layer(x=1:max_iter, y=model(1:max_iter), Geom.line, Theme(default_color=color("red"))),
		Guide.title("Energy vs Metropolis Iteration"),
		Guide.XLabel("Iteration"), Guide.YLabel("Energy")), "equilibration")

#=
 = Part d - Average energy and variation
=#
const conv_time = 4*equil_time

conv_energies = energies[conv_time+1:end]
conv_average = mean(conv_energies)
conv_variation = mean((conv_energies .- conv_average) .^ 2)
println("Average energy after convergence: ", conv_average)
println("Variation after convergence: ", conv_variation)

#=
 = Part e - correlation time
=#
const corr_time = 2500
const max_m = round(Int, (max_iter - conv_time)/2)

corrs = zeros(max_m)
for m = 1:max_m
	corrs[m] = mean((conv_energies[m+1:end] .- conv_average).*(conv_energies[1:end-m] .- conv_average))/conv_variation
end

model(t) = exp(-t/corr_time)

save_plot(plot(layer(x=1:max_m, y=corrs, Geom.line),
	layer(x=1:max_m, y = model(1:max_m), Geom.line, Theme(default_color=color("red"))),
	Guide.title("Correlation vs Separation"), Guide.XLabel("Separation (Iterations)"),
		Guide.YLabel("Correlation")), "correlation")

#=
 = Part f - Boltzmann distribution
=#
const num_groups = 1000

samples = []
append!(samples, particles)
for group in 1:num_groups
	for iter = 1:4*corr_time
		metropolis(particles)
	end
	append!(samples, particles)
end


save_plot(plot(layer(x=samples, Geom.histogram(density=true)), layer(x=0:0.01:5, y=2*e.^(-(0:0.01:5)*2), Geom.line,
	Theme(default_color=color("red"))), Guide.title("Particle Location PDF"), Guide.XLabel("Location"),
		Guide.YLabel("Probability Density")), "boltzmann")