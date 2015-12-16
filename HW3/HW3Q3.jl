#Including reusable code
include("HW3.jl")

#Constants of the problem
const alpha = 3/2
const r0 = 1
const Ns = 100
const Nc = 500
const L = 250
const N = 400
const d = L/N

function RL_walk(length)
  #Seeding at a random point in the box
  x = L*rand(3)
  pts = Vector{Float64}[]
  push!(pts,x)
  #Walking
  for step in 1:length-1
    #Get next distance
    r = r0*(1-rand())^(-1/alpha)
    #Get next direction
    theta, phi = 2*pi*rand(2)
    dir = [cos(theta)*sin(phi),sin(theta)*sin(phi),cos(phi)]
    #Adding the step
    x += r * dir
    x = mod(x,L) #Periodic boundary conditions
    push!(pts,x)
  end
  return pts
end

function CIC(particles)
  grid = zeros(N,N,N)
  for pt in particles
    p = pt/d #De-dimensionalizing
    for i in [1, -1]
      for j in [1,-1]
        for k in [1,-1]
          x = (floor(p[1] + (i+1)/2))
          y = (floor(p[2] + (j+1)/2))
          z = (floor(p[3] + (k+1)/2))
          dens = (1-abs(x-p[1]))*(1-abs(y-p[2]))*(1-abs(z-p[3]))
          grid[x%N+1,y%N+1,z%N+1] += dens
        end
      end
    end
  end
  return grid
end

function power_spectrum(F, m = N)
  PS = zeros(int(m/2))
  nums = zeros(int(m/2))
  for i = 1:m/2 + 1
    for j = 1:m/2 + 1
      for k = 1:m/2 + 1
        idx = floor(sqrt((i-1)^2 + (j-1)^2 + (k-1)^2))+1
        if idx > m/2
          break
        end
        PS[idx] += abs2(F[i,j,k])
        nums[idx] += 1
      end
    end
  end
  return PS ./ nums ./ m^3
end

function power_spectrum_corrected(F, m = N)
  PS = zeros(int(m/2))
  nums = zeros(int(m/2))
  for i = 1:m/2 + 1
    for j = 1:m/2 + 1
      for k = 1:m/2 + 1
        idx = floor(sqrt((i-1)^2 + (j-1)^2 + (k-1)^2))+1
        if idx > m/2
          break
        end
        w = abs2(sin(pi*i/N)*sin(pi*j/N)*sin(pi*k/N)/(pi*(i*j*k)/N))
        PS[idx] += abs2(F[i,j,k]) / w^2
        nums += 1
      end
    end
  end
  return PS ./ nums ./ m^3
end

#=
 = Part a - Rayleigh-Levy random walk
=#
particles = Vector{Float64}[]
for n in 1:Nc
  append!(particles, RL_walk(Ns))
end

#=
 = Part b - Particle density grid
=#
grid = CIC(particles)

#=
 = Part c - FFT
=#
gridF = rfft(grid)
println(sum(grid)) #Checking to make sure we still have 50,000 particles

#=
 = Part d - Power spectrum
=#
PS = power_spectrum(gridF)
save_plot(plot(x = 2:size(PS)[1], y = PS[2:end], Scale.x_log10, Scale.y_log10,
               Guide.title("Power Spectrum of Particle Density from RL walks"),
               Guide.XLabel("Wavenumber (Fundamentals)"), Guide.YLabel("Power (Density^2)")), "q3ps")

#Correcting
PS_c = power_spectrum_corrected(gridF)
save_plot(plot(x = 2:size(PS_c)[1], y = PS_c[2:end], Scale.x_log10, Scale.y_log10,
               Guide.title("Corrected Power Spectrum of Particle Density from RL walks"),
               Guide.XLabel("Wavenumber (Fundamentals)"), Guide.YLabel("Power (Density^2)")), "q3ps_c")
