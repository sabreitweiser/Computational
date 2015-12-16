#Including reusable code
include("HW3.jl")

const weeks = 104 #For use later

#Reading the data and putting into arrays
#Use the data frame reader to import data
data = readtable("TandH.dat", names = [:Temp, :Hum], separator = ' ', header = false)
Temp = array(data[:Temp]) #Temperature
Hum = array(data[:Hum]) #Humidity


#=
 = Part a - Means and Std Devs
 = I'm not saying it's trivial, but...
=#

Temp1 = Temp[1:weeks/2] #Year 1 Temp
Temp2 = Temp[weeks/2 + 1:end] #Year 2 Temp
Hum1 = Hum[1:weeks/2] #Year 1 Hum
Hum2 = Hum[weeks/2 + 1:end] #Year 2 Hum
println("Average Temp: ", mean(Temp))
println("Average Hum: ", mean(Hum))
println("Std Dev Temp: ", std(Temp))
println("Std Dev Hum: ", std(Hum))
println("Average Temp Year 1: ", mean(Temp1))
println("Average Hum Year 1: ", mean(Hum1))
println("Std Dev Temp Year 1: ", std(Temp1))
println("Std Dev Hum Year 1: ", std(Hum1))
println("Average Temp Year 2: ", mean(Temp2))
println("Average Hum Year 2: ", mean(Hum2))
println("Std Dev Temp Year 2: ", std(Temp2))
println("Std Dev Hum Year 2: ", std(Hum2))

#=
 = Part b - FFTs
=#

#Taking the FFT of temp
TempF = rfft(Temp)
#Drop first point because it's just the average
save_plot(plot(x = 1:weeks/2, y = abs2(TempF[2:end] ./ 104),
               Guide.title("Temperature Power Spectrum"),
               Guide.XLabel("Wavenumber (multiple of fundamental)"),
               Guide.YLabel("Power (Celsius^2)")), "temp_PS")

#Taking FFT of Humidity
HumF = rfft(Hum)
save_plot(plot(x = 1:weeks/2, y = abs2(HumF[2:end] ./ 104),
               Guide.title("Humidity Power Spectrum"),
               Guide.XLabel("Wavenumber (multiple of fundamental)"),
               Guide.YLabel("Power (%^2)")), "hum_PS")

#=
 = Part c - Fluctuation comparison
=#
#Julia indexes from 1, so 3 corresponds to the n = 2 term
#Set to 0 to remove power from seasonal changes (2 years in our period)
TempF[3] = 0
Temp_ns = irfft(TempF, weeks)
println("Temp std/mean, seasonal variation removed: ", std(Temp_ns)/mean(Temp_ns))
save_plot(plot(x = 1:weeks, y = Temp_ns, Guide.title("Temperature (Without Seasons)"),
               Guide.XLabel("Week"), Guide.YLabel("Temperature (Celsius)")), "temp_ns")

HumF[3] = 0
Hum_ns = irfft(HumF, weeks)
println("Hum std/mean, seasonal variation removed: ",std(Hum_ns)/mean(Hum_ns))
save_plot(plot(x = 1:weeks, y = Hum_ns, Guide.title("Humidity (Without Seasons)"),
               Guide.XLabel("Week"), Guide.YLabel("Humidity (%)")), "hum_ns")
