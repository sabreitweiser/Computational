#Importing reusable code and constants
include("HW2.jl")

#=
 = Part a - RK4 solver
=#

#Nothing to put here really, other than to refer to the RK4 implementation in HW2.jl

#=
 = Part b - Dimensionless Equations
=#
#Dimensional constants
const M0 = M
const r0 = 1e3
const e0 = 3.0006e-3 * M * c^2 / r0^3

#"Fixed" Equation of State
# Julia complains when I raise a negative number to a fractional power;
# this only became a problem when I increased p0 in part d, but I'm unsure why.
# I tried to correct by raising the absolute value and then fixing the sign,
# which makes sense if you think of it as raising by the numerator and then
# taking the denominator root;
# but I'm not sure if this is the right way to handle this.
# It seems to give the right answer, probably because small negative
# pressure is just a numerical error anyways
eps(p) = 2.4216*(abs(p)^(3/5))*abs(p)/p + 2.8663*p

#Mass diffeq
dm(r,p,m) = (r0^3*e0/(M0*c^2))*4*pi*r^2*eps(p)

#Newtonian differential equation
dp_nr(r,p,m) = -(G*M0/(c^2*r0))*eps(p)*m/r^2
#TOV differential equation
function dp_r(r,p,m)
  ret = -(G*M0/(c^2*r0))
  ret *= eps(p)/r^2
  ret *= (1+(p/eps(p)))
  ret *= (m+(e0*r0^3/(M0*c^2))*(4*pi*r^3*p))
  ret *= 1/(1-(G*M0/(c^2*r0))*(2*m/r))
  ret
end

#Wrapping the into combined diff eqs
dpm_nr(r, pm) = [dp_nr(r,pm[1],pm[2]),dm(r, pm[1], pm[2])]
dpm_r(r, pm) = [dp_r(r,pm[1],pm[2]),dm(r, pm[1], pm[2])]

#=
 = Part c - Newtonian vs Relativistic plotting
=#

#Helper function to return the orbit
#  (Inefficient, but best way to get orbital positions in an array)
function walk(f, start, step, initial; cutoff = 1e-7)
  R = Float64[start]
  r = start
  P = Float64[initial[1]]
  Mass = Float64[initial[2]]
  val = initial
  while val[1] > cutoff
    val = rk4(f, r, r + step, 1, val)
    r += step
    push!(R, r)
    push!(P, val[1])
    push!(Mass, val[2])
  end
  return R,P,Mass
end

start = 1e-7
step = 0.1
initial = [0.01, 0.0]
R_nr, P_nr, M_nr = walk(dpm_nr, start, step, initial)
R_r, P_r, M_r = walk(dpm_r, start, step, stop, initial)
save_plot(plot(layer(x=R_nr, y=P_nr, Geom.line), layer(x=R_r, y=P_r, Geom.point),
          Guide.title("Pressure vs Radius, Newtonian vs TOV"),
          Guide.XLabel("Radius (km)"), Guide.YLabel("Pressure (p0)")), "pvr")
save_plot(plot(layer(x=R_nr, y=M_nr, Geom.line), layer(x=R_r, y=M_r, Geom.point),
          Guide.title("Mass vs Radius, Newtonian vs TOV"),
          Guide.XLabel("Radius (km)"), Guide.YLabel("Mass (M0)")), "mvr")

#=
 = Part d - Maximum Mass
=#
#Obtaining parametric values for R and M
Rs = Float64[]
Ms = Float64[]
for p0 = 0.01:0.005:0.10
  initial = [p0, 0.0]
  R_r, P_r, M_r = walk(dpm_r, start, step, stop, initial)
  println(R_r[end]," ", M_r[end])
  push!(Rs, R_r[end])
  push!(Ms, M_r[end])
end

#Finding the maximum Mass and the corresponding Radius
max = indmax(Ms)
println(Ms[max], " solar masses")
println(Rs[max], " km")
#=prints
0.776381828189917 solar masses
10.30000009999998 km
=#

#Plotting
save_plot(plot(x=Rs, y=Ms,
               Guide.title("Mass vs Radius"), Guide.XLabel("Radius (km)"),
               Guide.YLabel("Mass (M0)")), "MvR2")
