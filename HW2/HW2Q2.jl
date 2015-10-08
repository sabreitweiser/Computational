#Importing reusable code and constants
include("HW2.jl")

#Constants
const delta = 0.05
const r0 = h^2 / (G*M)

#=
 = Part a - orbital equation
=#
du(t,u) = [u[2], (G*M/h^2)*(u[1]*r0)^delta - u[1]]

#=
 = Part b - computing the shift
=#
step = 2.0^-9
val = [1/perihelion, 0.0]

#First step it to 3*pi/2, so we've passed aphelion
#Fun fact, I only know what "aphelion" and "perhilion"
#  are because of kerbal space program
phi = 3*pi/2
val = rk4(du, 0, phi, phi/step, val)

#Now look for the derivative to become negative
# indicating we've passed perihelion
Phi = Float64[phi]
U = Float64[val[1]]
while true
  val = rk4(du, phi, phi + step, 1, val)
  phi += step
  push!(Phi, phi)
  push!(U, val[1])
  if val[2] < 0
    break
  end
end

#=
 = This part copied from Q1 with only small changes
=#
#Chopping off the last few points near perihilion
phi = Phi[end-2:end]
parab = U[end-2:end]

#Fitting to a parabola a*phi^2 + b*phi + x = u (standard matrix inversion)
coeffs = [(phi .^ 2).'; phi.'; ones(phi).'].'
A,B,C = inv(coeffs)*parab

#Parabola minimum at -b/(2*a), subtract off 2*pi to get shift
shift = -B/(2*A) - 2*pi

#Calculating shift in arcseconds per century
#  3600 arcseconds in degree
#  180/pi degrees in radian
#  365/87 orbits in year
#  100 years in century
println(shift, " radians per orbit")
println(shift*3600*365*100/87*180/pi, " arcseconds per century")
#=prints
0.1664540101656291 radians per orbit
1.4404328180334458e7 arcseconds per century
=#

#Plotting (for use in report)
phi2 = phi[1]:step/100:phi[end]
p(x) = A*x^2 + B*x + C
parab2 = map(p, phi2)
save_plot(plot(layer(x=phi, y = 1 ./ parab, Geom.point),
               layer(x=phi2, y = 1 ./ parab2, Geom.line),
               Guide.title("Force Law Perihelion"), Guide.XLabel("phi"),
               Guide.YLabel("r (m)")),"Question 2 Perihelion")
