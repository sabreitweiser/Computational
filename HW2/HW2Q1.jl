#Importing reusable code and constants
include("HW2.jl")

#=
 = Part a - RK4 algorithm
=#
#See the HW2.jl file for the RK4 implementation

#Diff eq, non-relativistic
du_nr(t, u) = [u[2], 1/(a*(1-e^2)) - u[1]]

#=
 = Part b - Non-relatistic RK4 vs Exact comparison
=#

#Helper function to return the orbit
#  (Inefficient, but best way to get orbital positions in an array)
function orbit(f, start, step, stop, initial)
  Phi = Float64[start]
  phi = start
  Ret = Float64[initial[1]]
  val = initial
  while phi < stop
    val = rk4(f, phi, phi + step, 1, val)
    phi += step
    push!(Phi, phi)
    push!(Ret, val[1])
  end
  return Phi, Ret
end

#Checking for closed orbit
Phi, numerical = orbit(du_nr, 0, 0.1, 4*pi, [1/perihelion, 0.0])
save_plot(plot(x = cos(Phi) ./ numerical, y = sin(Phi) ./ numerical,
          Guide.title("RK4 orbit, step = 0.1"), Guide.XLabel("x (m)"),
               Guide.YLabel("y (m)")), "rk4orbit")

#Exact Newtonian eq
u(phi) = (1+e*cos(phi))/(a*(1-e^2))

#Function to orbit and plot
function non_rel_orbit(step)
  name = string("Runge Kutta vs Exact, step = ", step, " (non-relativistic)")
  save_name = string("rk4Vexact", integer(100*step))
  exact = Float64[]
  Phi, numerical = orbit(du_nr, 0, step, 4*pi, [1 / perihelion, 0])

  for phi in Phi
    push!(exact, u(phi))
  end

  uvphi = plot(layer(x=Phi, y=exact, Geom.line),
               layer(x=Phi, y=numerical, Geom.point),
               Guide.title(name), Guide.XLabel("phi"),
               Guide.YLabel("u = 1/r (1/m)"))
  save_plot(uvphi, save_name)
end

#Using a few different step sizes
non_rel_orbit(0.1)
non_rel_orbit(0.2)
non_rel_orbit(0.05)

#=
 = Part c - Calculating the perihelion shift
=#

#Diff eq, relativistic
du_r(t,u) = [u[2], G*M/h^2 + 3*G*M/(c^2)*u[1]^2 - u[1]]

#Step through the orbit using our relatavistic diff eq
step = 2.0^-10
pt = 2*pi
Phi, calc = orbit(du_r, 0, step, pt+step, [1 / perihelion, 0])

#Chopping off the last few points near 2*pi
phi = Phi[end-2:end]
parab = calc[end-2:end]

#Fitting to a parabola a*phi^2 + b*phi + x = u (standard matrix inversion)
coeffs = [(phi .^ 2).'; phi.'; ones(phi).'].'
A,B,C = inv(coeffs)*parab

#Parabola minimum at -b/(2*a), subtract off 2*pi to get shift
shift = -B/(2*A) - 2*pi

println(shift, " radians per orbit")
#Calculating shift in arcseconds per century
#  3600 arcseconds in degree
#  180/pi degrees in radian
#  365/87 orbits in year
#  100 years in century
println(shift*3600*365*100/87*180/pi, " arcseconds per century")
#=prints
5.075586999581105e-7 radians per orbit
43.92229467890694 arcseconds per century
=#

#Plotting (for use in report)
phi2 = phi[1]:step/100:phi[end]
p(x) = A*x^2 + B*x + C
parab2 = map(p, phi2)
save_plot(plot(layer(x=phi, y = 1 ./ parab, Geom.point),
               layer(x=phi2, y = 1 ./ parab2, Geom.line),
               Guide.title("Relativistic Perihelion"), Guide.XLabel("phi"),
               Guide.YLabel("r (m)")),"Relativistic Perihelion")
