"""

"""
function infected(h, he, vstatus, istatus)
    vertices = getvertices(h, he)
    sum = 0
    for v in vertices
        if istatus[v.first] == 0
            sum += vstatus[v.first]
        end
    end
    sum
end


"""

"""
function f(ninfected, c)
    if ninfected > c
        return c
    else
        return ninfected
    end
end


"""

"""
function TVHSIS(
        df,
        intervals,
        node_index_map,
        he_index_map,
        Δ,
        δ;
        mytitle = "",
        vstatus = nothing,
        per_infected = 20,
        c=12,
        τd=0.2,
        τₕ=0.2,
        τᵥ=0.20,
        γₕ=0.05,
        γᵥ=0.05,
        αᵥ = 0.05,
        βᵥ = 0.0,
        αₑ = 0.05,
        βₑ = 0.0,
        ntrials=1,
        immunizenode = false,
        immunizeedge = false
    )

    per_infected_sim = Dict{Int, Array{Float64, 1}}()
    per_infected_sim_trials = Dict{Int, Array{Float64, 1}}()

    if c == nothing
        fname = "AAMAS/results/1year-Δ:$(Δ)-perinfected:$(per_infected)-c:median-τₕ:$(τₕ)-τᵥ:$(τᵥ)-γₕ:$(γₕ)-γᵥ:$(γᵥ)-αᵥ:$(αᵥ)-βᵥ:$(βᵥ).csv"
    else
        #fname = "AAMAS/results/tokyo/hg/1year-Δ:$(Δ)-perinfected:$(per_infected)-c:$(c)-τₕ:$(τₕ)-τᵥ:$(τᵥ)-γₕ:$(γₕ)-γᵥ:$(γᵥ)-αᵥ:$(αᵥ)-βᵥ:$(βᵥ)-αᵥ:$(αₑ)-βᵥ:$(βₑ).csv"
        fname = "src/experiments/results/53/csv/$(mytitle).csv"
    end

    #hgs = Array{Any, 1}()

    open(fname, "a") do fhandle
        write(
            fhandle,
            "time,delta,newusers,movedusers,c,td,th,tv,gh,gv,aq,bq,aeq,beq,density,avg_degree,perc_tot_infected,tot_he_infected,tot_he_infected_1\n")

        for trial=1:ntrials
            println(trial)

            #random status initialization
            h = nothing
            added, moved = 0, 0

            usersepoc = rand(0:0, 1, length(keys(node_index_map)))

            if vstatus == nothing
                vstatus = rand(1:1, 1, length(keys(node_index_map)))
                vrand = rand(0:100, 1, length(keys(node_index_map)))
                for i=1:length(keys(node_index_map))
                    if per_infected  <= vrand[i]
                        vstatus[i] = 0
                    end
                end
            end
            hstatus = rand(0:0, 1, length(keys(he_index_map)))

            #immunization
            istatus = rand(0:0, 1, length(keys(node_index_map)))
            irand = rand(0:100, 1, length(keys(node_index_map)))

            ihestatus = rand(0:0, 1, length(keys(he_index_map)))
            iherand = rand(0:100, 1, length(keys(he_index_map)))

            vnextstatus = copy(vstatus)
            nextistatus = copy(istatus)

            henextstatus = copy(hstatus)
            nextihestatus = copy(ihestatus)

            # 1 infetto
            # vstatus = rand(0:0, 1, length(keys(node_index_map)))
            # vstatus[rand(1:length(keys(node_index_map)))] = 1

            push!(
                get!(per_infected_sim, 0, Array{Float64, 1}()),
                sum(vstatus)/length(vstatus))

            push!(
                get!(per_infected_sim_trials, 0, Array{Float64, 1}()),
                sum(vstatus)/length(vstatus))

            #
            #
            #simulate
            #
            #
            #println("0 -- Vstatus $(vstatus)")
            #println("istatus $(istatus)")
            #println("0 -- Hstatus $(hstatus)")
            for t=1:length(intervals)
                h, a, m = generatehg(h, df, get(intervals, t, 0).first, get(intervals, t, 0).second, node_index_map, he_index_map, t, true, usersepoc)

                #h = buildhg(df, get(intervals, t, 0).first, get(intervals, t, 0).second, convert(Dates.Millisecond, Dates.Hour(Δ)), node_index_map, he_index_map)

                #push!(hgs, copy(h))
                #println(vstatus)

                if t == 200
                    if immunizenode
                        for i=1:length(keys(node_index_map))
                            if αᵥ > irand[i]/100
                                nextistatus[i] = 1
                            end
                        end
                    elseif immunizeedge
                        for i=1:length(keys(he_index_map))
                            if αₑ > iherand[i]/100
                                nextihestatus[i] = 1
                            end
                        end
                    end
                end

                if h == nothing
                    continue
                end

                #push!(hgs, copy(h))
                dist = Array{Int, 1}()
                if c == nothing
                    for he=1:size(h)[2]
                        push!(dist, length(getvertices(h, he)))
                    end
                    c = median(dist)
                    println(t, " -- ", c)
                end

                #number of possibile interactions in the hypergraph
                density = 0
                for he=1:size(h)[2]
                    density += length(getvertices(h, he))
                end

                #hg avg node degree
                avgdegree = 0
                for v=1:nhv(h)
                    avgdegree += length(gethyperedges(h, v))
                end
                avgdegree /= length(node_index_map)

                hstatus_greater_1 = 0
                for he=1:size(h)[2]
                    if hstatus[he] == 1 && length(getvertices(h, he)) > 1
                        hstatus_greater_1+=1
                    end
                end

                for iter=1:1
                    #
                    # direct influence
                    #
                    avgdirectdegree = 0
                    for v=1:size(h)[1]
                        if usersepoc[v] == 1

                                i = 0
                                for he in gethyperedges(h, v)
                                    for u in getvertices(h, he.first)
                                        if v != u.first
                                            if abs(h[v, he.first] - h[u.first, he.first]) <= δ.value
                                                #println("$(h[v, he.first]) -- $(h[u.first, he.first]) -- $(δ)")
                                                if vstatus[v] == 0
                                                    i += vstatus[u.first]
                                                end
                                                avgdirectdegree += 1
                                            end
                                        end
                                    end
                                end
                                if vstatus[v] == 0 && rand(1)[1] < 1 - ℯ ^ -(τd * i)
                                    vnextstatus[v] = 1
                                end

                        end
                    end

                    #println("$(t) -- AvgDirectDegree $(avgdirectdegree) / $(sum(usersepoc)) = $(avgdirectdegree \= sum(usersepoc))")

                    avgdirectdegree \= sum(usersepoc)



                    #
                    # undirect influence
                    #
                    for he=1:size(h)[2]
                        if ihestatus[he] == 1
                            continue
                        end
                        if length(getvertices(h, he)) > 1 && hstatus[he] == 0
                            i = infected(h, he, vstatus, istatus)
                            if rand(1)[1] <  1 - ℯ ^ -(τₕ * f(i, c))
                                #hstatus[he] = 1
                                henextstatus[he] = 1
                            #    if t > 100 && (rand(1)[1] < 1 - ℯ ^ - αₑ)
                                    #ihestatus[he] = 1
                            #        nextihestatus[he] = 1
                            #    end
                            end
                        elseif rand(1)[1] <  1 - ℯ ^ - γₕ
                            #hstatus[he] = 0
                            henextstatus[he] = 0
                            #ihestatus[he] = 0
                        # elseif ihestatus[he] == 1 && rand(1)[1] < 1 - ℯ ^ - (βₑ) #out of quarentine
                        #     ihestatus[he] = 0
                        #elseif t > 100 && (rand(1)[1] < 1 - ℯ ^ - αₑ)
                            #ihestatus[he] = 1
                        #    nextihestatus[he] = 1
                        end
                    end

                    for v=1:size(h)[1]
                        if istatus[v] == 1 #se sono immunizzato
                            continue
                        end
                        if usersepoc[v] == 1
                            if vstatus[v] == 0
                                i = 0
                                for he in gethyperedges(h, v)
                                    if length(getvertices(h, he.first)) > 1
                                        if ihestatus[he.first] == 0
                                            i += hstatus[he.first]
                                        end
                                    end
                                end
                                if rand(1)[1] < 1 - ℯ ^ -(τᵥ * f(i, c))
                                    #vstatus[v] = 1
                                    vnextstatus[v] = 1
                                    #go to quarentine with prob p
                        #            if t > 100 && (rand(1)[1] < 1 - ℯ ^ - αᵥ)
                                        #istatus[v] = 1
                                        #vstatus[v] = 0
                        #                nextistatus[v] = 1
                        #                vnextstatus[v] = 0
                        #            end
                                end
                            elseif rand(1)[1] < 1 - ℯ ^ - γᵥ #healthy and out of quarentine
                                    #vstatus[v] = 0
                                    vnextstatus[v] = 0
                                    #istatus[v] = 0
                            # elseif istatus[v] == 1 && rand(1)[1] < 1 - ℯ ^ - (βᵥ) #out of quarentine
                            #         istatus[v] = 0
                        #    elseif t > 100 && (rand(1)[1] < 1 - ℯ ^ - αᵥ)
                                #istatus[v] = 1
                                #vstatus[v] = 0
                        #        nextistatus[v] = 1
                        #        vnextstatus[v] = 0
                            end
                        end
                    end
                    vstatus = copy(vnextstatus)
                    istatus = copy(nextistatus)

                    hstatus = copy(henextstatus)
                    ihestatus = copy(nextihestatus)

                    # push!(
                    #     get!(per_infected_sim, iter+(t-1)*100, Array{Float64, 1}()),
                    #     sum(vstatus)/(length(vstatus) - sum(istatus))
                    # )

                    push!(
                        get!(per_infected_sim, iter+(t-1)*1, Array{Float64, 1}()),
                        sum(vstatus)/(length(vstatus) - sum(istatus))
                    )

                    #println("$(t) -- Vstatus $(vstatus)")
                    #println("$(t) -- istatus $(istatus)")
                    #println("$(t) -- Hstatus $(hstatus)")

                    write(fhandle,"$(t),$(Δ),$(a),$(m),$(c),$(τₕ),$(τᵥ),$(γₕ),$(γᵥ),$(αᵥ),$(βᵥ),$(αₑ),$(βₑ),$(density),$(avgdegree),$(sum(vstatus)/length(vstatus)),$(sum(hstatus)),$(hstatus_greater_1)\n")
                end
            end
        end
    end
    per_infected_sim#, hgs
end
