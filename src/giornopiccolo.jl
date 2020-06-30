using TVHEpidemicDynamics
using DataFrames
using CSV
using Dates
using PyPlot
using Statistics




dati="C:\\Users\\Utente\\Desktop\\tesi uni\\TVHEpidemicDynamics.jl\\src\\weeplace_checkins.csv" ##fonte

header=[:userid, :venueid, :UTCtime,  :lat, :lng, :city, :category]
dateformat = "yyyy-mm-ddTHH:MM:SS"


df = CSV.read(
    dati;
    copycols = true,
    header = header,
    dateformat = dateformat
    )

mindate = minimum(df[!, :UTCtime])
print(mindate)

##2003-11-02T06:04:47 E' IL PIù PICCOLO è il primo. 