#=
 = Including reusable code
=#
include("HW5.jl")

#=
 = Setup
=#
# Constant parameters of the problem
const L = 50
const size = L^2

#n is usually less than 1 million, so cut it off by 10 million
const max_n = 10000000
n = 0

#number of samples in part b
const num_sample = 5

#Helper function for periodic BCs
idx(i) = mod(i-1, size) + 1

#Calculate beta from temperature, using de-dimensionalized units (k=J=1)
beta_calc(T) = 1/T

# IMPLEMENTATION NOTE: Spins are stored as bit arrays for speed and memory efficiency

#Calculate magnetization from bit array
magnetization(spins) = abs(size - 2*sum(spins))

#hot seed(random)
seedhot() = bitrand(size)

#cold seed (all spin down or up, chosen randomly)
function seedcold()
	if bitrand()[1]
		return trues(size)
	else
		return falses(size)
	end
end

#=
 = Part a - Metropolis
=#
function metropolis(spins, beta)
	pos = floor(Int,rand()*size)
	spin = spins[idx(pos)]
	dE = -4*((spin $ spins[idx(pos-1)]) + (spin $ spins[idx(pos+1)])
				+ (spin $ spins[idx(pos-L)]) + (spin $ spins[idx(pos+L)])) + 8
	#dE = +-2, so -4*(spin1 xor spin2) + 2 for each
	if dE <= 0
		spins[idx(pos)] = ~spin
	else
		if rand() < exp(-beta*dE)
			spins[idx(pos)] = ~spin
		end
	end
end

T_test = 2.5 #Temperature to test at, very close to critical
beta_test = beta_calc(T_test)

spins_hot = seedhot()
spins_cold = seedcold()

mag_hot = Float64[]
mag_cold = Float64[]

#Find the value at which they converge (are within the expected sqrt(size) variation)
for step = 1:max_n
	metropolis(spins_hot, beta_test)
	metropolis(spins_cold, beta_test)
	push!(mag_hot, magnetization(spins_hot))
	push!(mag_cold, magnetization(spins_cold))
	if abs(magnetization(spins_hot)-magnetization(spins_cold)) < sqrt(L)
		n = step
		break
	end
end
save_plot(plot(layer(x=1:n, y=mag_hot, Geom.line), layer(x=1:n, y=mag_cold, Geom.line),
		Guide.title("Hot vs Cold Initialization Magnetization, T=2.5"), Guide.XLabel("Iteration"),
		Guide.YLabel("Abs(Magnetization)")), "hotcold")
println(n, " Steps needed for convergence at T=", T_test)

#Now assign a value for the number of iterations needed for convergence
const num_iter = 1000000


#=
 = Part b - Magnetization
=#

function converge(spins, beta)
	for step = 1:num_iter
		metropolis(spins, beta)
	end
end

const num_converge = 1
const num_group = 1

function predict_mag(T)
	println(T)
	sum = 0.0
	beta = beta_calc(T)
	for i = 1:num_sample
		spins = seedhot()
		for j = 1:num_converge
			converge(spins, beta)
		end
		for j = 1:num_group
			converge(spins, beta)
			sum += magnetization(spins)
		end
	end
	return sum/(num_sample*num_group)
end

model(T) = int(T .< 2.269)

temps = 0:0.05:3.0
mags = map(predict_mag, temps)
save_plot(plot(layer(x=temps, y=mags/size, Geom.point),
	layer(x=temps, y = model(temps), Geom.line, Theme(default_color=color("red"))),
	Guide.title("Abs(Magnetization) vs Temperature, Hot Initialization"),
	Guide.XLabel("Temperature (J/k)"), Guide.YLabel("Abs(Magnetization)/Size")), "tempvsmag_hot")