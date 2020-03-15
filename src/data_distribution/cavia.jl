#cerco i chekins solo negli spazi diretti, cosndiero contatti diretti.

using CSV
using DataFrames
using Dates
using SimpleHypergraphs

export evaluateintervals!, evaluatedensity!
export evaluatedistance!, evaluate_direct_contacts_distribution!
export evaluate_location_distribution!

export buildparamsdirect, buildhg, generatehg
export TVHSIS

include("helper.jl")
include("epidemics/TVH.jl")
include("epidemics/TVHSIS.jl")
using DataFrames
using CSV
using Dates
using PyPlot
using Statistics

#atom mi  dà problemi, ho risolto in questo semplice modo
dataset = "C:/Users/Utente/Documents/GitHub/TVHEpidemicDynamics.jl/src/data_distribution/Brightkite_totalCheckins.txt"
header = [:userid, :UTCtime, :lat, :lng,:venueid]
#dateformat = "e u d H:M:S +0000 Y" #y-m-dTH:M:SZ
dateformat= "y-m-dTH:M:SZ"
df = CSV.read(
    dataset;
    copycols = true,
    header = header,
    dateformat = dateformat
    )


df=dropmissing(df)



#getting first checkin within the data, alcuni non hanno data e quindi skippo perchè mi va in errore.
mindate = minimum((df[!, :UTCtime]))
minmonday = Dates.DateTime("2009-05-12T00:00:00")
mintime = minmonday

#
# DISTRIBUTION PARAMETERS
#
Δₕ = 24*7
Δₘ = convert(Dates.Millisecond, Dates.Hour(Δₕ))

# Evaluating the total number of time slots
# within the observation period (last_checkin_date - minmonday)
# using Δₘ as discretization param
nintervals = ceil(Int,(maximum(skipmissing(df[!, :UTCtime])) - minmonday)/Δₘ)

# interval_id -> (start_date, end_date)
intervals = Dict{Int, Pair{DateTime, DateTime}}()
evaluateintervals!(mintime, maximum(skipmissing(df[!, :UTCtime])), intervals, nintervals, Δₘ)


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
δₕ = 4
δₘ = convert(Dates.Millisecond, Dates.Hour(δₕ))

# Evaluating the total number of time slots
# within the week presenting the majority of checkins
# using δₘ as discretization parameter
nintervals_δ = ceil(Int,(maxweek.second - maxweek.first)/δₘ)

# interval_id -> (start_date, end_date)
intervals_δ = Dict{Int, Pair{DateTime, DateTime}}()
evaluateintervals!(maxweek.first, maxweek.second, intervals_δ, nintervals_δ, δₘ)

# Evaluatin checkin distribution within each interval.
distancewithinintervals = Dict{Int, Array{Any, 1}}()
evaluatedistance!(intervals_δ, distancewithinintervals, df)



"""
    Plot checkin distribution within 7th-13th May.
"""
clf()

fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(111)

ylabel("Time in Seconds", fontweight="semibold")
xlabel("Time Intervals", fontweight="semibold")

for t=1:length(keys(distancewithinintervals))
    values = get(distancewithinintervals, t, Array{Any, 1}())
    ax.boxplot(values, widths=0.5, positions=[t])
end

ax.set_xlim(0, 43)

xlabels = ["12 May", "04", "08", "12", "16", "20", "13 May", "04", "08", "12", "16", "20", "14 May", "04", "08", "12", "16", "20", "15 May", "04", "08", "12", "16", "20", "16 May", "04", "08", "12", "16", "20", "17 May", "04", "08", "12", "16", "20", "18 May", "04", "08", "12", "16", "20"]
#xlabels = ["7\nMay", "", "", "12h", "", "", "8\nMay", "", "", "12h", "", "", "9\nMay", "", "", "", "", "", "10\nMay", "", "", "", "", "", "11\nMay", "", "", "", "", "", "12\nMay", "", "", "", "", "", "13\nMay", "", "", "", "", ""]

xticks(1:42, xlabels, rotation=80)

gcf()

plt.tight_layout(.5)

savefig("src/data_distribution/plots/checkin_distributiontestGAETANO.png")



"""
    Plot checkin distribution within 7th-13th May.
"""
toplot = Vector{Array{Any,1}}()
cut = 0

for t=1:length(keys(distancewithinintervals))
    values = get(distancewithinintervals, t, Array{Any, 1}())
    #println(t, "-->", values)
    #cut = 60
    # if length(values) == 0 || t==16
    #     cut = 60
    # else
    #     cut = quantile(values, 0.9)
    # end

    # try
    #     println("$(t) -- $(cut) -- $(median(values)) -- $(length(values[values.<cut])/length(values)) -- $(length(values[values.>cut])/length(values))")
    # catch ArgumentError
    #     println("$(t) -- 60 -- NaN -- $(length(values[values.<cut])/length(values)) -- $(length(values[values.>cut])/length(values))")
    # end

    #filter!(x -> x < cut, values)
    push!(toplot, values)
end

clf()

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(111)
ax.boxplot(toplot)

xlabels = "12 May", "04", "08", "12", "16", "20", "13 May", "04", "08", "12", "16", "20", "14 May", "04", "08", "12", "16", "20", "15 May", "04", "08", "12", "16", "20", "16 May", "04", "08", "12", "16", "20", "17 May", "04", "08", "12", "16", "20", "18 May", "04", "08", "12", "16", "20"]
xticks(1:42, xlabels, rotation=80)

ylabel("Time in Seconds")
xlabel("Time Intervals")

gcf()

plt.tight_layout(.5)

savefig("src/data_distribution/plots/directcontacts_distributionGAETANO.png")
