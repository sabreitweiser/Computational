#=
 = Including reusable code
=#
include("HW4.jl")

#=
 = Part a - Galerkin Method
=#
# Constant Parameters
const M = 20
const t_min = 0
const ht = 0.1
const t_max = 200
t_grid = [0 , 20, 40, 80, 120, 200]
t_len = size(t_grid)[1]
const alpha = 1.0
const L = 40

# Helper functions
k(n) = pi*n/(L)
idx(n) = n < 0 ? 2*M+2+n : n+1


f = zeros(Complex{Float64}, 2*M+1)
f[idx(1)] = 1/16
f[idx(-1)] = 1/16

function deriv(f)
	ret = zeros(f)
	for n = -M:M
		ret[idx(n)] += im*alpha*k(n)^3*f[idx(n)]
		for np = max(n-M,-M):min(M,n+M)
			ret[idx(n)] -= f[idx(n-np)]*f[idx(np)]*im*k(np)
		end
	end
	return ret
end #deriv

function val(f,x)
	ret = 0.0
	for n = -M:M
		ret += real(f[idx(n)]*exp(im*k(n)*x))
	end
	return ret
end #val

d(t, f) = deriv(f)
t_last = t_min
coeffs = zeros(Complex{Float64}, t_len, size(f)[1])
PS = Vector{Float64}[]
for i = 1:t_len
	t = t_grid[i]
	coeffs[i, :] = f
	push!(PS, abs2(f[1:M]))
	f = rk4(d, t_last, t, (t-t_last)/ht, f)
end

#=
 = Part b - Plots
=#
save_plot(plot([(x) -> val(coeffs[i, :], x) for i = 1:t_len], -L, L, color=t_grid,
	Guide.title("Solution Evolving"), Guide.YLabel("f"), Guide.XLabel("x"),
	Guide.colorkey("Time")),"q4_solution")

#=
 = Part c - Power Spectra
=#

save_plot(plot(x=1:M, y=PS[1], Guide.title("Power Spectrum t=0"), Guide.YLabel("Power"),
	Guide.XLabel("Wavenumber"), Geom.line),"q4_ps_0")
save_plot(plot(x=1:M, y=PS[2], Guide.title("Power Spectrum t=20"), Guide.YLabel("Power"),
	Guide.XLabel("Wavenumber"), Geom.line),"q4_ps_20")
save_plot(plot(x=1:M, y=PS[3], Guide.title("Power Spectrum t=40"), Guide.YLabel("Power"),
	Guide.XLabel("Wavenumber"), Geom.line),"q4_ps_40")
save_plot(plot(x=1:M, y=PS[4], Guide.title("Power Spectrum t=80"), Guide.YLabel("Power"),
	Guide.XLabel("Wavenumber"), Geom.line),"q4_ps_80")
save_plot(plot(x=1:M, y=PS[5], Guide.title("Power Spectrum t=120"), Guide.YLabel("Power"),
	Guide.XLabel("Wavenumber"), Geom.line),"q4_ps_120")
save_plot(plot(x=1:M, y=PS[6], Guide.title("Power Spectrum t=200"), Guide.YLabel("Power"),
	Guide.XLabel("Wavenumber"), Geom.line),"q4_ps_200")

#=
 = Part d - Explanation
=#
# For the alpha = 0.1 graphs, I just lowered alpha manually and saved the graphs,
#	to save myself some coding. If you want to reproduce this you can easily change
#	alpha in the constants block at the top of this file.
