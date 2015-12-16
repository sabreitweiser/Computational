#Including reusable code
include("HW3.jl")

const l = 100 #Grid size

#Reuseable code to compute power spectrum fron 2d FT
function power_spectrum(F, m = l)
  PS = zeros(int(m/2))
  nums = zeros(int(m/2))
  for i = 1:m/2 + 1
    for j = 1:m/2 + 1
      idx = floor(sqrt((i-1)^2 + (j-1)^2))+1
      if idx > m/2
        continue #Drop points in shells that go beyond the domain
      end
      PS[idx] += abs2(F[i,j])
      nums[idx] += 1
    end
  end
  return PS ./ nums / m^2
end

for n = 1:3 #For each of the grids

  #Reading the data in
  filename = string("gaussian", n, ".dat")
  data = readtable(filename, separator = ' ', names = [:x, :y, :f], header = false)
  f_grid = reshape(array(data[:f]), l, l)

#=
 = Part a - FFT
=#
  F = rfft(f_grid)

#=
 = Part b - power spectrum
=#
  PS = power_spectrum(F)
  save_plot(plot(x = 2:size(PS)[1], y = PS[2:end], Scale.x_log10, Scale.y_log10,
                 Guide.Title(string("Power Spectrum Number ", n)),
                 Guide.XLabel("Wavenumber"), Guide.YLabel("Power")), string("ps", n))

#=
 = Part c - power spectrums
=#
  #No code needed here

#=
 = Part d - resampling
=#
  f_grid = f_grid[1:2:l, 1:2:l]
  F = rfft(f_grid)
  PS = power_spectrum(F, l/2)
  save_plot(plot(x = 2:size(PS)[1], y = PS[2:end], Scale.x_log10, Scale.y_log10,
                 Guide.Title(string("Resampled Power Spectrum Number ", n)),
                 Guide.XLabel("Wavenumber"), Guide.YLabel("Power")), string("ps2", n))
end
