"""

"""


"""

"""
function buildparamsdirect(
    fname::AbstractString,
    columns::AbstractVector,
    dateformat::AbstractString;
    Δ::Millisecond = Dates.Millisecond(14400000),
    δ::Millisecond = Dates.Millisecond(1000),
    limit::Int = typemax(Int),
    datelimit::Union{DateTime, Nothing} = nothing,
    mindatelimit::Union{DateTime, Nothing} = nothing
    )

    df = CSV.read(
        fname;
        copycols=true,
        header=columns,
        dateformat=dateformat,
        limit=limit)

    # df = CSV.read(fname;
    #     copycols=true, header=[:userid, :UTCtime, :lat, :lng, :venueid],
    #     dateformat="y-m-dTH:M:SZ", #e u d H:M:S +0000 Y
    #     limit=limit)

    df = dropmissing(df, :UTCtime)

    #println(size(df))

    if datelimit != nothing
        df = filter(r -> r.UTCtime <= datelimit, df)
    end

    if mindatelimit != nothing
        df = filter(r -> r.UTCtime >= mindatelimit, df)
    end

    sort!(df, :UTCtime)

    #vertices
    users = unique(df, :userid)
    nusers = nrow(users)

    #hyperedges - upper buond
    places = unique(df, :venueid)
    nplaces = nrow(places)

    #time slots
    nintervals = ceil(Int,(maximum(df[!, :UTCtime]) - minimum(df[!, :UTCtime]))/Δ)
    intervals = Dict{Int, Pair{DateTime, DateTime}}()

    #evaluate timeframes
    mintime = minimum(df[!, :UTCtime])
    for i=1:nintervals
        offset = mintime + Δ + δ > maximum(df[!, :UTCtime]) ?  maximum(df[!, :UTCtime]) + Dates.Millisecond(1) :  mintime + Δ + δ
        #get!(intervals, i, Pair{DateTime, DateTime}(convert(Dates.DateTime, (mintime-delta)), convert(Dates.DateTime, (offset+delta))))
        get!(intervals, i, Pair{DateTime, DateTime}(convert(Dates.DateTime, (mintime - δ)), convert(Dates.DateTime, (offset))))

        mintime = offset - δ
    end

    #mapping node -> user
    node_index_map = Dict{String, Int}()
    index = 1
    #add vertices
    for user in eachrow(users)
        userid = string(user[:userid])
        if !haskey(node_index_map, userid)
            push!(node_index_map, userid => index)
            index += 1
        end
    end

    #mapping he -> place
    he_index_map = Dict{String, Int}()
    index = 1
    #add he
    for place in eachrow(places)
        placeid = string(place[:venueid])
        if !haskey(he_index_map, placeid)
            push!(he_index_map, placeid => index)
            index += 1
        end
    end

    #println(node_index_map)
    #println(he_index_map)

    df, intervals, node_index_map, he_index_map
end


"""

"""
function buildhg(
    df,
    mindate,
    maxdate,
    delta,
    node_index_map,
    he_index_map
    )

    prevdf = filter(r -> ((r.UTCtime >= mindate-delta) && (r.UTCtime < (mindate))), df)
    currdf = filter(r -> ((r.UTCtime >= (mindate)) && (r.UTCtime < (maxdate))), df)
    nextdf = filter(r -> ((r.UTCtime >= (maxdate))) && (r.UTCtime < maxdate+delta), df)

    if size(currdf)[1] == 0
        return nothing
    end

    if size(prevdf)[1] != 0
        toadd_p = enrichdf(currdf, prevdf, delta)
        for r in toadd_p
            push!(currdf, r)
        end
    end
    if size(nextdf)[1] != 0
        toadd_n = enrichdf(currdf, nextdf, delta)
        for r in toadd_n
            push!(currdf, r)
        end
    end

    h = Hypergraph{Int, String, String}(length(keys(node_index_map)), length(keys(he_index_map)))

    for user in node_index_map
        set_vertex_meta!(h, user.first, user.second)
        #println(user.first, "=>", user.second)
    end

    for place in he_index_map
        set_hyperedge_meta!(h, place.first, place.second)
        #println(place.first, "=>", place.second)
    end

    for checkin in eachrow(currdf)
        # se un utente visita lo stesso luogo nello stesso timeframe?
        # viene memorizzato il timestamp dell'ultima volta che lo visita
        setindex!(h, Dates.value(checkin.UTCtime), get(node_index_map, string(checkin.userid), -1), get(he_index_map, checkin.venueid, -1))
    end

    h
end


"""

"""
function enrichdf(df1, df2, delta)
    toreturn = Array{Any, 1}()

    for r1 in eachrow(df2)
        for r2 in eachrow(df1)
            if (r1.venueid == r2.venueid) && (abs(r1.UTCtime - r2.UTCtime) < delta)
                push!(toreturn, r1)
            end
        end
    end

    toreturn
end


"""
"""
function generatehg(h, df, mindate, maxdate, node_index_map, he_index_map, t, all, usersepoc)
    added = 0
    moved = 0

    vadded = rand(0:0, 1, length(keys(node_index_map)))

    #initialize hg
    if h == nothing
        h = Hypergraph{Int, String, String}(length(keys(node_index_map)), length(keys(he_index_map)))

        for user in node_index_map
            set_vertex_meta!(h, user.first, user.second)
            #println(user.first, "=>", user.second)
        end

        for place in he_index_map
            set_hyperedge_meta!(h, place.first, place.second)
            #println(place.first, "=>", place.second)
        end
    else
        for v_index=1:nhv(h)
            max = 0
            id = 0
            for he in gethyperedges(h, v_index)
                if getindex(h, v_index, he.first) > max
                    max = getindex(h, v_index, he.first)
                    id = he.first
                end
                setindex!(h, nothing, v_index, he.first)
            end
            if id != 0
                setindex!(h, max, v_index, id)
            end
        end
    end

    currdf = filter(r -> ((r.UTCtime >= (mindate)) && (r.UTCtime < (maxdate))), df)

    for checkin in eachrow(currdf)
        #if the user was not in any other venue, ie hyperedge, it is a new node
        #hg_userid = get(node_index_map, string(checkin.userid), -1))
        if usersepoc[get(node_index_map, string(checkin.userid), -1)] == 0
            added += 1
            usersepoc[get(node_index_map, string(checkin.userid), -1)] = 1
        else
            moved += 1
        end

        # se un utente visita lo stesso luogo nello stesso timeframe?
        # viene memorizzato il timestamp dell'ultima volta che lo visita

        if t == 1 || all
            # try
            setindex!(h, Dates.value(checkin.UTCtime), get(node_index_map, string(checkin.userid), -1), get(he_index_map, string(checkin.venueid), -1))
            # catch e
            #     println(checkin.userid, "-->", get(node_index_map, string(checkin.userid), -1), "    ",
            #         checkin.venueid, "-->", get(he_index_map, checkin.venueid, -1), "    ", size(h))
            #     showerror(stdout, e)
            #     println()
            # end
        end
    end
    h, added, moved

end
