"""
Creates a hypergraph from a dataset of checkins.
"""


"""
Loads a dataframe from a file `fname` containing checkin-data.
For each checkin stores user and place in two dictionaries.
Creates a dictionary of time intervals describing the evolution of the relationships
between users and places.
"""
function buildparamsdirect(fname::AbstractString,
    columns::AbstractVector,
    dateformat::AbstractString;
    Δ::Millisecond = Dates.Millisecond(14400000),
    δ::Millisecond = Dates.Millisecond(1000),
    limit::Int = typemax(Int),
    mindate::Union{DateTime,Nothing} = nothing,
    maxdate::Union{DateTime,Nothing} = nothing)

    df = CSV.read(fname;
    copycols = true,
    header = columns,
    dateformat = dateformat,
    limit = limit)

    df = dropmissing(df, :UTCtime)
    df = filter(r->(isnothing(mindate) || r.UTCtime >= mindate) &&
                    (isnothing(maxdate) || r.UTCtime <= maxdate), df)
    sort!(df, :UTCtime)

    # vertices
    users = unique(df, :userid)
    numusers = nrow(users)

    # hyperedges - upper buond
    places = unique(df, :venueid)
    numplaces = nrow(places)

    # time slots
    mintime = minimum(df[!, :UTCtime])
    maxtime = maximum(df[!, :UTCtime])

    nintervals = ceil(Int, (maxtime - mintime) / Δ)
    intervals = Dict{Int,Pair{DateTime,DateTime}}()

    # evaluate timeframes
    for i = 1:nintervals
        offset = mintime + Δ + δ > maxtime ? maxtime + Dates.Millisecond(1) :
        mintime + Δ + δ
        get!(intervals, i, Pair{DateTime,DateTime}(convert(Dates.DateTime,
        (mintime - δ)), convert(Dates.DateTime, (offset))))

        mintime = offset - δ
    end

    # mapping node -> user
    node_index_map = Dict{String,Int}()
    index = 1

    # add vertices
    for user in eachrow(users)
        userid = string(user[:userid])
        if !haskey(node_index_map, userid)
            push!(node_index_map, userid => index)
            index += 1
        end
    end

    # mapping he -> place
    he_index_map = Dict{String,Int}()
    index = 1
    # add he
    for place in eachrow(places)
        placeid = string(place[:venueid])
        if !haskey(he_index_map, placeid)
            push!(he_index_map, placeid => index)
            index += 1
        end
    end

    df, intervals, node_index_map, he_index_map
end


"""
Returns a hypergraph for a given set of vertices and hyperedges based on the time
parameters.
"""
function buildhg(df,
    mindate,
    maxdate,
    δ,
    node_index_map,
    he_index_map)

    prevdf = filter(r->((r.UTCtime >= mindate - δ) && (r.UTCtime < mindate)), df)
    currdf = filter(r->((r.UTCtime >= mindate) && (r.UTCtime < maxdate)), df)
    nextdf = filter(r->((r.UTCtime >= maxdate)) && (r.UTCtime < maxdate + δ), df)

    if size(currdf)[1] == 0
        return nothing
    end

    enrichdf!(prevdf, currdf, δ)
    enrichdf!(nextdf, currdf, δ)

    h = inithg(currdf, node_index_map, he_index_map)
    h
end


"""
Adds every checkin `c1` in `df1` to `df2` if any checkin `c2` in `df2` has been
made in the same place of `c1` and in a close period of time
"""
function enrichdf!(df1, df2, δ)
    toadd = Array{Any,1}()

    for r1 in eachrow(df1)
        for r2 in eachrow(df2)
            if (r1.venueid == r2.venueid) && (abs(r1.UTCtime - r2.UTCtime) < δ)
                push!(toadd, r1)
            end
        end
    end

    for row in toadd
        push!(df2, row)
    end
end


"""
Returns a hypergraph for a given set of vertices and hyperedges specifying the
number of new nodes
"""
function generatehg(h,
    df,
    mindate,
    maxdate,
    node_index_map,
    he_index_map,
    t,
    all,
    usersepoc)

    added = 0
    moved = 0
    vadded = rand(0:0, 1, length(keys(node_index_map)))

    currdf = filter(r->((r.UTCtime >= mindate) && (r.UTCtime < maxdate)), df)

    # initialize hg
    if isnothing(h)
        h = inithg(currdf, node_index_map, he_index_map)
    else
        for v_index = 1:nhv(h)
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

    for checkin in eachrow(currdf)
        # if the user was not in any other venue, ie hyperedge, it is a new node
        usersepoc[get(node_index_map, string(checkin.userid), -1)] == 0 ?
        usersepoc[get(node_index_map, string(checkin.userid), -1)] = 1 : moved += 1

        # if a user visits the same place in the same timeframe
        # the timestamp of the last time he visits it is stored
        if t == 1 || all
            setindex!(h, Dates.value(checkin.UTCtime), get(node_index_map,
            string(checkin.userid), -1), get(he_index_map, string(checkin.venueid), -1))
        end
    end
    added = size(currdf)[1] - moved

    h, added, moved
end

"""
Creates a hypergraph from a dataframe of checkins for a given set of vertices and
hyperedges
"""
function inithg(df, node_index_map, he_index_map)
    h = Hypergraph{Int,String,String}(length(keys(node_index_map)), length(keys(he_index_map)))

    for user in node_index_map
        set_vertex_meta!(h, user.first, user.second)
    end

    for place in he_index_map
        set_hyperedge_meta!(h, place.first, place.second)
    end

    for checkin in eachrow(df)
        # if a user visits the same place in the same timeframe
        # the timestamp of the last time he visits it is stored
        setindex!(h, Dates.value(checkin.UTCtime), get(node_index_map,
        string(checkin.userid), -1), get(he_index_map, checkin.venueid, -1))
    end

    h
end
