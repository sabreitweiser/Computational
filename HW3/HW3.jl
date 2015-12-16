#=
 = Reusable code file
=#
using Gadfly
using DataFrames

cd(dirname(@__FILE__)) #Interactive interpreter is weird about files, setting the directory maually

#Function to save plots for use in report
function save_plot(plot, name)
  folder_path = "C:\\Users\\Alex\\Documents\\Computational\\HW3\\"
  draw(PNG(string(folder_path, name, ".png"),6inch, 4inch), plot)
end