#Import our plotting library
using Gadfly

#Function to plot integral error as a function of step
function i_plot(X, Y, title)
  return plot(x=X, y=Y, Guide.title(title), Guide.XLabel("Number of intervals"),
              Guide.YLabel("Integral Rel. Error"), Scale.x_log10, Scale.y_log10)
end #function i_plot

#Function to save plots for use in report
function save_plot(plot, name)
  folder_path = "C:\\Users\\Alex\\Documents\\Computational\\HW1\\"
  draw(PNG(string(folder_path, name, ".png"),6inch, 4inch), plot)
end #save_plot

#Function returning the relative error, all using single precision
function rel_err(act::Float32, pred::Float32)
  return abs((act-pred)/act)
end #rel_err

#Constant for two in single precision
const TWO_SP = float32(2)

#Function that implements the trapezoid rule for integration
function trap_integral(f, N::Unsigned, x_min::Float32, x_max::Float32)
  h = (x_max - x_min)/(N-1)
  tot = (f(x_min)+f(x_max))/TWO_SP
  for i = 1:N-2
    tot += f(x_min + i*h)
  end
  tot *= h
  return tot
end #trap_integral

#Constants for three, four in single precision
const THREE_SP = float32(3)
const FOUR_SP = float32(4)

#Function that implements the Simpson rule for integration
function simp_integral(f, N::Unsigned, x_min::Float32, x_max::Float32)
  h = (x_max - x_min)/(N-1)
  tot = (f(x_max) + f(x_min))
  for i = 1:2:N-2
      tot += FOUR_SP * f(x_min + i*h)
  end
  for i = 2:2:N-3
      tot += TWO_SP * f(x_min + i*h)
  end
  tot *= h / THREE_SP
  return tot
end #simp_integral

#Function that implements the Gaussian-Legendre Quadrature
function gauss_integral(f, N::Unsigned, x_min::Float32, x_max::Float32)
  alpha = (x_max - x_min)/TWO_SP
  beta = (x_max + x_min)/TWO_SP
  X, W = Base.gauss(Float32, N)
  tot = float32(0)
  for (x,w) in zip(X,W)
    tot += w*f(alpha*x + beta)
  end
  tot *= alpha
end #gauss_integral

#Setting up our function and its analytical integral
f(x::Float32) = exp(-x)
actual = float32(1-exp(-1))

#Helper lines
calc_n(int_f, n) = rel_err(actual, int_f(f, n, float32(0), float32(1)))
N = unsigned(2.^(1:20))#NB: period indicates vectorized operation, like matlab

#Calculating and plotting Trapezoid rule integration relative errors, for N = powers of two + 1
calc_trap_np1(n) = calc_n(trap_integral, n+1)
trap_error = map(calc_trap_np1, N)
save_plot(i_plot(N, trap_error, "Trapezoid integration error"), "trap")

#Calculating and plotting Simpson rule integration relative errors, for N = powers of two + 1
calc_simp_np1(n) = calc_n(simp_integral, n+1)
simp_error = map(calc_simp_np1, N)
save_plot(i_plot(N, simp_error, "Simpson integration error"), "simp")

#Calculating and plotting Gauss-Lagrange quadrature integration relative error, for N = 1 to 100
N = 1:100
calc_gauss_n(n) = calc_n(gauss_integral, n)
gauss_error = map(calc_gauss_n, N)
save_plot(i_plot(N, gauss_error, "Gaussian Quadrature integration error"), "gauss")
