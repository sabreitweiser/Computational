#=
 = Including reusable code
=#
include("HW4.jl")

#=
 = Setup
=#
# Constant parameters of the problem
const v_max = 95.33 #speed limit, ft/s (~65mph)
const rho_max = 5.0 #max density, 1 car every 20ft
const L = 52800.0 #10 mile track, in ft
const T = 3600 #One hour, in seconds
const hx = 100.0 #ten foot spacing, half a car length
const ht = 0.1 #One tenth of a second timing
const x_grid = -L/2:hx:L/2
const x_len = size(x_grid)[1]
const t_grid = 0:ht:T
const t_len = size(t_grid)[1]

#=
 = Part a - Intuition
=#
# Defining relevant functions
v(rho) = v_max*(1-rho/rho_max)
c(rho) = v_max*(1-2*rho/rho_max)
f(rho) = rho*v(rho)

#=
 = Part b - Initial Conditions
=#
# Setting up the initial conditions
rho = zeros(t_len,x_len)
rho[1,round(Int64,x_len/4):round(Int64,x_len/2)] = rho_max
#plot(y = rho[1,:], x=x_grid/5280, Geom.line)

#=
 = Part c - Solving
=#
const mu = ht/hx
# FTCS
function ftcs_step(rho)
	ret = zeros(rho)
	for x in 1:x_len
		xp1 = mod(x,x_len)+1
		xm1 = mod(x-2,x_len)+1
		ret[x] = rho[x] - mu/2*(f(rho[xp1]) - f(rho[xm1]))
	end
	return ret
end #ftcs
for step in 1:t_len-1
	rho[step+1, :] = ftcs_step(rho[step, :])
end
println(sum(rho[t_len,:]))
for time = 0:5
	t = round(Integer, (t_len-1)*time / 5)
	title = string("Traffic after ", round(Integer, t*ht/60), " min (FTCS)")
	save_plot(plot(layer(y = rho[t+1,:], x=x_grid/5280, Geom.line),
		Guide.title(title), Guide.XLabel("Location, miles"),
		Guide.YLabel("Car Density, Cars per 100 ft"),
		layer(y=zeros(rho[1,:]), x=x_grid/5280, Geom.line, Theme(default_color=color("grey")))),
		title)
end


# Resetting for Lax
rho = zeros(t_len,x_len)
rho[1,round(Int64,x_len/4):round(Int64,x_len/2)] = rho_max*0.99

# Lax
function lax_step(rho)
	ret = zeros(rho)
	for x in 1:x_len
		xp1 = mod(x,x_len)+1
		xm1 = mod(x-2,x_len)+1
		ret[x] = 1/2*(rho[xp1] + rho[xm1]) - mu/2*(f(rho[xp1]) - f(rho[xm1]))
	end
	return ret
end #lax
for step in 1:t_len-1
	rho[step+1, :] = lax_step(rho[step, :])
end
println(sum(rho[t_len,:]))
for time = 0:5
	t = round(Integer, (t_len-1)*time / 5)
	title = string("Traffic after ", round(Integer, t*ht/60), " min (Lax)")
	save_plot(plot(layer(y = rho[t+1,:], x=x_grid/5280, Geom.line),
		Guide.title(title), Guide.XLabel("Location, miles"),
		Guide.YLabel("Car Density, Cars per 100 ft"),
		layer(y=zeros(rho[1,:]), x=x_grid/5280, Geom.line, Theme(default_color=color("grey")))),
		title)
end