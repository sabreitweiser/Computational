#Import our plotting library
using Gadfly

#Function to plot properly
function my_plot(X, Y, title, ylabel)
  return plot(x=X, y=Y, Guide.title(title), Guide.XLabel("Number of steps"), Guide.YLabel(ylabel))
end #function my_plot

#Function to save plots
function save_plot(plot, name)
  folder_path = "C:\\Users\\Alex\\Documents\\Computational\\HW1\\"
  draw(PNG(string(folder_path, name, ".png"),6inch, 4inch), plot)
end

#Defining constants of the LCG
const A = 9301
const C = 49297
const M = 233280

#Random number generator generator, using LCG algo
function random_gen(seed::Integer)
  state = seed
  function random()
    state = (A*state + C) % M
    return state / float32(M) #Need to promote operation to floating point
  end #random
  return random
end #random_gen

#Random walk function; iteratively produces position from +-1 random walk
function random_walk(seed::Integer, n_max::Integer)
  next = random_gen(seed) #Get an RNG
  x = 0 #Start the walk at 0
  function walk()
    for step = 1:n_max
      if next() > 0.5 #50/50 chance to walk forward or back
        x += 1
      else
        x -= 1
      end
      produce(x) #yield the current position
    end
  end
  @task walk() #pre-wrap the producer in a Task
end #random_walk

squared = zeros(Int64, 500)
cubed = zeros(Int64, 500)
fourth = zeros(Int64, 500)
const WALK_LENGTH = 500
seeds = 1:1000
for seed in seeds
  i = 1
  for x in random_walk(seed, WALK_LENGTH)
    squared[i] += x^2
    cubed[i] += x^3
    fourth[i] += x^4
    i += 1
  end
end

#Correcting
squared /= 1000

cubed /= 1000
cubed ./= squared.^(3.0/2.0)

fourth ./= (squared).^2
fourth /= 1000
fourth -= 3

N = 1:WALK_LENGTH

#Plotting
sq = my_plot(N, squared, "sigma^2", "sigma^2")
tr = my_plot(N, cubed, "s3", "s3")
fr = my_plot(N, fourth, "s4", "s4")

#Saving
save_plot(sq, "sq")
save_plot(tr, "tr")
save_plot(fr, "fr")
