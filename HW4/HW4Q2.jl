#=
 = Including reusable code
=#
include("HW4.jl")

#=
 = Some setup
=#
#Setting up constant parameters of the problem
const x_max = 1.5
const x_min = 0
const y_max = 2
const y_min = 0
const grid_size = 100
const hx = (x_max - x_min)/grid_size
const hy = (y_max - y_min)/grid_size
const x_grid = linspace(x_min,x_max,grid_size)
const y_grid = linspace(y_min,y_max,grid_size)

LOWER_BC(x) = 1.5 - 2*x
UPPER_BC(x) = 2.75 - 1.5*x
#=
save_plot(plot([LOWER_BC, UPPER_BC],0,1.5, Guide.yticks(ticks=[y_min:grid_size:y_max]),
	Guide.xticks(ticks=[x_min:grid_size:x_max]), Guide.title("Region of interest"),
	Guide.XLabel("x"), Guide.YLabel("y")), "Q2_region")
=#

#=
 = Part a - GS Relaxation
=#

#Constants used in Gauss-Seidel
const hx2 = hx^2
const hy2 = hy^2
const h1 = 1/(2*(hx2+hy2))
const h2 = hx2*hy2
const prec = 10
num_iter = convert(Integer,prec*grid_size^2/4)


phi = zeros(grid_size, grid_size)
function relax(phi)
	max_diff = 0
	for x = 2:grid_size-1
		xval = x_grid[x]
		y_bottom = max(2, convert(Integer,ceil(LOWER_BC(xval)/hy)))
		y_top = min(grid_size-1, convert(Integer,floor(UPPER_BC(xval)/hy)))
		for y = y_bottom:y_top
			last = phi[x,y]
			phi[x,y] = h1*(hy2*(phi[x+1,y] + phi[x-1,y]) + hx2*(phi[x,y+1] + phi[x,y-1])) + (h1*h2)
			diff = abs(last - phi[x,y])
			if diff > max_diff
				max_diff = diff
			end
		end
	end
	return max_diff
end #relax

#=
 = Part b - Error (and actual relaxation)
=#
Err = Float64[]
b = 0
for iter = 1:num_iter
	err = relax(phi)
	if err == 0
		b = iter-1
		break
	end
	push!(Err, err)
end
save_plot(plot(x=1:b, y=Err, Scale.y_log10, Geom.line,
	Guide.title("Error vs Number of Iterations (GS)"), Guide.XLabel("Iterations"),
	Guide.YLabel("Error")), "q2_gs_error")

#=
 = Part c - Contour Plot
=#
save_plot(plot(z=phi, x=x_grid, y=y_grid, Geom.contour,
	Guide.title("Solution (GS)"), Guide.colorkey("Phi")),
	"q2_gs_solution")

#=
 = Part d - SOR Method
=#
const w = 1.2
phi = zeros(grid_size, grid_size)
function relax(phi)
	max_diff = 0
	for x = 2:grid_size-1
		xval = x_grid[x]
		y_bottom = max(2, convert(Integer,ceil(LOWER_BC(xval)/hy)))
		y_top = min(grid_size-1, convert(Integer,floor(UPPER_BC(xval)/hy)))
		for y = y_bottom:y_top
			last = phi[x,y]
			phi[x,y] = (1-w)*phi[x,y] + w*h1*(hy2*(phi[x+1,y] + phi[x-1,y]) + hx2*(phi[x,y+1] + phi[x,y-1])) + (h1*h2)
			diff = abs(last - phi[x,y])
			if diff > max_diff
				max_diff = diff
			end
		end
	end
	return max_diff
end #relax

Err = Float64[]
for iter = 1:b
	err = relax(phi)
	push!(Err, err)
end

save_plot(plot(x=1:b, y=Err, Scale.y_log10, Geom.line,
	Guide.title("Error vs Number of Iterations (SOR, w=1.2)"), Guide.XLabel("Iterations"),
	Guide.YLabel("Error")), "q2_sor_error")


save_plot(plot(z=phi, x=x_grid, y=y_grid, Geom.contour,
	Guide.title("Solution (SOR, w=1.2)"), Guide.colorkey("Phi")),
	"q2_sor_solution")