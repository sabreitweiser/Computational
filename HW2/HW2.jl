#=
 = Reusable code file
=#
using Gadfly

#Function to save plots for use in report
function save_plot(plot, name)
  folder_path = "C:\\Users\\Alex\\Documents\\Computational\\HW2\\"
  draw(PNG(string(folder_path, name, ".png"),6inch, 4inch), plot)
end

#Global constants, in MKS
const G = 6.674e-11 #m^3/(kg*s^2)
const M = 1.989e30 #kg
const c = 2.998e8 #m/s
const e = 0.2056 #unitless
const a = 57.91e9 #m
const perihelion = a*(1-e)
const h = 58.98e3 * a*(1-e) #Max velocity times perihelion distance (constant), m^2/s

#4th order explicit Runge-Kutta algorithm
function rk4(f, t0, tN, N, u0)
  h = (tN-t0)/N
  u = u0
  t = t0
  for i = 1:N
    k1 = f(t,u)
    k2 = f(t+(h/2),u+(h/2)*k1)
    k3 = f(t+(h/2),u+(h/2)*k2)
    k4 = f(t+h, u+h*k3)
    u .+= (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    t += h
  end
  return u
end #rk4