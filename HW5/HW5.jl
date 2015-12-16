#=
 = Reusable code file
=#
using Gadfly

#Function to save plots for use in report
function save_plot(plot, name)
  folder_path = "C:\\Users\\Alex\\Documents\\Computational\\HW5\\"
  draw(PNG(string(folder_path, name, ".png"),6inch, 4inch), plot)
end