#2003-11-02T06:04:47
using TVHEpidemicDynamics
using DataFrames
using CSV
using Dates
using PyPlot
using Statistics

"""
 Studying check-in weekly distribution
 within the same place, considering only direct contact.
"""

dataset = "C:\\Users\\Utente\\Desktop\\tesi uni\\TVHEpidemicDynamics.jl\\src\\weeplace_checkins.csv"
header=[:userid, :venueid, :UTCtime,  :lat, :lng, :city, :category]
dateformat = "yyyy-mm-ddTHH:MM:SS" #y-m-dTH:M:SZ

df = CSV.read(
    dataset;
    copycols = true,
    header = header,
    dateformat = dateformat
    )

#getting first checkin within the data
mindate = minimum(df[!, :UTCtime])

#evaluating data starting from the first monday
#minMonday = 2012-04-09T18:17:18
#minmonday = Dates.DateTime("2012-04-09T18:17:18")
minmonday = Dates.DateTime("2003-11-03T00:00:00")
mintime = minmonday


#
# DISTRIBUTION PARAMETERS
#
Δₕ = 24*30
Δₘ = convert(Dates.Millisecond, Dates.Hour(Δₕ))

# Evaluating the total number of time slots
# within the observation period (last_checkin_date - minmonday)
# using Δₘ as discretization param
print("aaaaaaaaaaaaaaaaaaAAAAAAAAAAAAAAAA\n")
print(typeof(df[!,:UTCtime]))
nintervals = ceil(Int,(maximum(df[!, :UTCtime]) - minmonday)/Δₘ)

# interval_id -> (start_date, end_date)
intervals = Dict{Int, Pair{DateTime, DateTime}}()
evaluateintervals!(mintime, maximum(df[!, :UTCtime]), intervals, nintervals, Δₘ)


# Count checkin data according to the discrete time intervals
# evaluated at the previous step
checkinsperinterval = Dict{Int, Int}()
evaluatedensity!(intervals, checkinsperinterval, df)

# Sort time intervals according to the number of checkins
# within each interval
sortedcheckinsperinterval = sort(collect(checkinsperinterval), by=x->x[2], rev=true)

# Working with the week containing
# the majority of checkins
maxweek = get(intervals, sortedcheckinsperinterval[1].first, nothing)
#maxweek = intervals[5]

print(maxweek)