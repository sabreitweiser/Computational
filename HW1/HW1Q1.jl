#Import our plotting library
using Gadfly

#Function to plot deriative error as a function of step
function d_plot(X, Y, title)
  return plot(x=X, y=Y, Guide.title(title), Guide.XLabel("Step"),
              Guide.YLabel("Derivative Rel. Error"), Scale.x_log10, Scale.y_log10)
end #function d_plot

#Function to save plots for use in report
function save_plot(plot, name, pname)
  folder_path = "C:\\Users\\Alex\\Documents\\Computational\\HW1\\"
  draw(PNG(string(folder_path, name, " ", pname, ".png"),6inch, 4inch), plot)
end

#Function for single precision derivative using forward difference algorithm
function forward_diff_derivative(f, x::Float32, h::Float32)
  return (f(x+h) - f(x))/h
end

#Julia interprets float literals as double precision, so creating 32-bit "two"
const TWO_SP = float32(2.0)

#Function for single precision derivative using center difference algorithm
#  Need to use TWO_SP so return expression doesn't get promoted to 64 bit
function center_diff_derivative(f, x::Float32, h::Float32)
  return (f(x + (h/TWO_SP)) - f(x - (h/(TWO_SP))))/h
end

#Julia interprets float literals as double precision, so creating sp constants
const THREE_SP = float32(3.0)
const FOUR_SP = float32(4.0)

#Function for single precision derivatice using extrapolated difference algorithm
#  Need to use single precision constants so the expression does not get promoted
function extrapol_diff_derivatice(f, x::Float32, h::Float32)
  return (FOUR_SP*center_diff_derivative(f, x, h/TWO_SP)
            - center_diff_derivative(f, x, h))/THREE_SP;
end

#Julia has a built in function to return the machine error!
const machine_error = eps(Float32)

#Function to compute relative error in single precision
function rel_err(predicted::Float32, actual::Float32)
  return abs((actual-predicted)/actual)
end

function plot_d_err(f, df, x::Float32, fname::String)
  pname = string(fname," ")
  if x < 1
    pname = string(pname, "01") #Don't want to confuse latex with a decimal before the picture extensions
  end
  if x > 1
    pname = string(pname, "10")
  end
  println(pname)
  #Actual value of derivative at that point (used for relative error)
  act = float32(df(x))

  #Pre-starting value for step, set to 16 to start at 8
  h = float32(16)
  #Array to track steps, used for plotting later
  H = Float32[]
  #Arrays to store relative erros
  f_err = Float32[]
  c_err = Float32[]
  e_err = Float32[]
  println("*** Calculating derivatives ***")
  #Decrease h by a factor of ten each time, keep going until h is less than the machine error
  while h > machine_error
    h *= float32(0.5)
    push!(H, h)

    #Calculating forward difference derivative and its relative error
    for_d = forward_diff_derivative(f, x, h)
    for_d_err = rel_err(for_d, act)
    push!(f_err, for_d_err)
    println("Forward ", "Step: ", h, " Error: ", for_d_err)

    #Calculating center difference derivative and its relative error
    cent_d = center_diff_derivative(f, x, h)
    cent_d_err = rel_err(cent_d, act)
    push!(c_err, cent_d_err)
    println("Center ", "Step: ", h, " Error: ", cent_d_err)

    #Calculating extraploated difference derivative and its relative errror
    ext_d = extrapol_diff_derivatice(f, x, h)
    ext_d_err = rel_err(ext_d, act)
    push!(e_err, ext_d_err)
    println("Extrap ", "Step: ", h, " Error: ", ext_d_err)
  end #while h > machine_error

  #Plotting
  println("*** Plotting ***")
  forward = d_plot(H, f_err, string("Forward Difference ", pname))
  center = d_plot(H, c_err,string("Center Difference ", pname))
  extrap = d_plot(H, e_err, string("Extrapolated Difference ", pname))

  #Printing Errors
  println("***Truncation Error slopes***")
  println("Forward: ", (log(f_err[5]) - log(f_err[4]))/(log(H[5]) - log(H[4])))
  println("Center: ", (log(c_err[5]) - log(c_err[1]))/(log(H[5]) - log(H[1])))
  println("Extrap: ", (log(e_err[5]) - log(e_err[1]))/(log(H[5]) - log(H[1])))

  println("***Roundoff Error slopes***")
  println("Forward: ", (log(f_err[end-5]) - log(f_err[end-1]))/(log(H[end-5]) - log(H[end-1])))
  println("Center: ", (log(c_err[end-5]) - log(c_err[end-1]))/(log(H[end-5]) - log(H[end-1])))
  println("Extrap: ", (log(e_err[end-5]) - log(e_err[end-1]))/(log(H[end-5]) - log(H[end-1])))

  save_plot(forward, "forward", pname)
  save_plot(center, "center", pname)
  save_plot(extrap, "extrap", pname)
end #function plot_d_err

df(x) = -sin(x)
plot_d_err(cos, df, float32(0.1), "cos")
plot_d_err(cos, df, float32(10), "cos")
plot_d_err(exp, exp, float32(0.1), "exp")
plot_d_err(exp, exp, float32(10), "exp")
