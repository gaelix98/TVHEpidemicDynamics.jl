using TVHEpidemicDynamics
using Dates
using CSV
using PyPlot

"""
    Experiment of section 5.3 of the paper
"""

#model params
Δₕ = 4
Δₘ = convert(Dates.Millisecond, Dates.Hour(Δₕ))

δₛ = 60#60*10
δₘ = convert(Dates.Millisecond, Dates.Minute(δₛ))

#lastcheckindate = Dates.DateTime("2013-4-03T18:17:18")
#firstcheckindate = Dates.DateTime("2012-05-07T18:17:18")
firstcheckindate = Dates.DateTime("2012-05-07T00:00:00")
lastcheckindate = Dates.DateTime("2012-06-07T00:00:00")

#firstcheckindate = Dates.DateTime("2009-09-04T05:17:38")

#lastcheckindate = Dates.DateTime("2013-4-03T18:17:18")

#Tokyo - datemin - 2013-4-03T18:17:18

#brightkite - datemax 2010-10-18T18:39:58
#brightkite - datemin 2008-03-21T20:36:21

#gowalla - datemax 2010-10-23T05:22:06
#gowalla - datemin 2009-02-04T05:17:38

#columns = [:userid, :UTCtime, :lat, :lng, :venueid]

#fname = "data/loc-gowalla_totalCheckins.txt"
fparams = "src/experiments/configs/53/aamas53.csv"

dataset = "data/dataset_TSMC2014_TKY.txt"
header = [:userid, :venueid, :catid, :catname, :lat, :lng, :timezoneoffset, :UTCtime]
dateformat = "e u d H:M:S +0000 Y" #y-m-dTH:M:SZ


c = 5

df, intervals, node_index_map, he_index_map =
    buildparamsdirect(
        dataset, header, dateformat;
        Δ=Δₘ, δ=δₘ,
        datelimit = lastcheckindate, mindatelimit=firstcheckindate)#, mindatelimit=firstcheckindate, limit=300000


paramsdf = CSV.read(
        fparams;
        copycols=true,
        header=[:data, :label, :per_infected, :Δ, :δ ,:βd, :βᵢ, :βₑ, :γₑ, :γₐ],#, :αᵥ, :βᵥ, :αₑ, :βₑ, :immunizenode, :immunizeedge],
        datarow=2)

test_data = Dict{String, Array{Any, 1}}()
for params in eachrow(paramsdf)
    push!(
        get!(test_data, SubString(params[:data], 1, 1), Array{Any, 1}()),
        params
    )
end


#mytitle = "aamasInfected60secs"
infected_distribution = Vector{Array{Float64,1}}()

# intervalsTest = Dict{Int, Pair{DateTime, DateTime}}()
# push!(
#     intervalsTest,
#     1 => Pair(Dates.DateTime("2012-05-07T17:17:00"), Dates.DateTime("2012-05-08T00:37:00"))
# )
# push!(
#     intervalsTest,
#     1 => Pair(Dates.DateTime("2012-05-07T21:17:00"), Dates.DateTime("2012-05-08T04:37:00"))
# )

#########################
#evaluate data
########################
intervals = unique(paramsdf, [:Δ, :δ])[!,[:Δ, :δ]]
intervals_data = Dict{String, Dict{Symbol, Any}}()

for i in eachrow(intervals)
    df, intervals, node_index_map, he_index_map =
        buildparamsdirect(
            dataset, header, dateformat;
            Δ=convert(Dates.Millisecond, Dates.Hour(i.Δ)), δ=convert(Dates.Millisecond, Dates.Minute(i.δ)),
            datelimit = lastcheckindate, mindatelimit=firstcheckindate)

    push!(
        get!(intervals_data, "$(i.Δ)$(i.δ)", Dict{Symbol, Any}()),
        :df => df,
        :intervals => intervals,
        :node_index_map => node_index_map,
        :he_index_map => he_index_map,
    )
end

###########
per_infected = 20

vstatus = rand(1:1, 1, size(unique(df, :userid))[1])
vrand = rand(0:100, 1, size(unique(df, :userid))[1])

for i=1:length(keys(node_index_map))
    if per_infected  <= vrand[i]
        vstatus[i] = 0
    end
end


###############################

infected_distribution = Vector{Array{Float64,1}}()

println(":data, :label, :per_infected, :Δ, :δ ,:βd, :βᵢ, :βₑ, :γₑ, :γₐ")

linestyles = ["solid", "dashed", "dashdot", "dotted"]
markers = ["", "", "", "", "x", "+"]

for testtype in keys(test_data)
    linestyle = 1
    marker = 1
    clf()
    figure(figsize=(7,4))
    labels = Array{String, 1}()
    mytitle = ""
    for test in get(test_data, testtype, nothing)

        println("$(test[:data]) $(test[:label]) $(test[:per_infected]) $(test[:Δ]) $(test[:δ]) $(test[:βd]) $(test[:βᵢ]) $(test[:βₑ]) $(test[:γₑ]) $(test[:γₐ])")

        runningparams = get(intervals_data, "$(test[:Δ])$(test[:δ])", Dict{Symbol, Any}())

        SIS_per_infected_sim =
            TVHSIS(
                get!(runningparams, :df, nothing),
                get!(runningparams, :intervals, nothing),
                get!(runningparams, :node_index_map, nothing),
                get!(runningparams, :he_index_map, nothing),
                test[:Δ],
                convert(Dates.Millisecond, Dates.Minute(test[:δ]));
                vstatus=vstatus,
                per_infected = test[:per_infected], c=c, τd=test[:βd],
                τₕ=test[:βᵢ], τᵥ=test[:βₑ], γₕ=test[:γₑ], γᵥ=test[:γₐ],
                #αᵥ=test[:αᵥ], βᵥ=test[:βᵥ], αₑ=test[:αₑ], βₑ=test[:βₑ],
                ntrials=1,
                mytitle=test[:data])

        #println(SIS_per_infected_sim)
        infected_distribution = Vector{Array{Float64,1}}()

        for t=1:length(keys(SIS_per_infected_sim))-1
            #println("$(t) -- $(get(SIS_per_infected_sim, t, Array{Float64, 1}()))")
            push!(infected_distribution, get(SIS_per_infected_sim, t, Array{Float64, 1}()))
        end

        #println(infected_distribution)
        #println("$(linestyles[linestyle]) --- $(markers[marker])")

        ylim(bottom=0.0, top=0.6)
        plot(infected_distribution, linestyle=linestyles[linestyle], marker=markers[marker], markevery=10, markersize=6.5)

        #title(gettitle(test))
        #ylabel("Δ = $(test[:Δ]) hours \n Infected nodes in %", labelpad=10, fontweight="semibold", fontsize="xx-large")
        ylabel("Infected nodes in %", fontweight="semibold", fontsize="x-large", labelpad=10)

        xlabel("Time intervals", fontweight="semibold", labelpad=10, fontsize="x-large")
        #title("δ = $(test[:δ]) minutes", pad=10, fontweight="semibold", fontsize="xx-large")

        tick_params(labelsize="large")


        #ax2 = ax.twiny()
        #ax2.set_xlabel("δ = 1 minute")

        varname = test[:data]#startswith(test[:data], "t") ? SubString(test[:data], 2) : test[:data]

        #push!(labels, "$(varname)=$(test[Symbol(varname)])")
        push!(labels, test[:label])
        mytitle = test[:data] #gettitle(test)
        #title(gettitle(test))
        #println(gettitle(test), " ------- ", varname, " ", test[Symbol(varname)])
        #println(gettitle(test), " ------- ", varname, " ", test[:data])
        #println(linestyle)
        linestyle = (linestyle + 1) % (length(linestyles)+1)
        marker = (marker + 1) % (length(markers)+1)

        if linestyle == 0
            linestyle = 1
        end
        if marker == 0
            marker = 1
        end
    end
    println(labels)
    legend(labels, fontsize="large", ncol=2)
    plt.tight_layout(.5)
    #savefig("src/AAMAS/resultplot/aamas/53/$(mytitle).png")
    savefig("src/experiments/results/53/plot/$(mytitle).png")
end

gcf()
