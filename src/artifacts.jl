using Artifacts

function _load_gr(filename)
    open(filename, "r") do file
        line = readline(file)
        first_str, second_str, first_num, second_num = split(line)
        n = parse(Int64, first_num)
        graph = SimpleGraph(n)

        for line in eachline(file)
            v1, v2 = split(line)
            add_edge!(graph, parse(Int64, v1), parse(Int64, v2))
        end
        return graph
    end
end

function graph_from_artifact(num::Int)
    try
        snum = lpad(string(num), 3, '0')
        file_name = joinpath(artifact"PACE2019", "vc-exact_" * snum * ".gr")
        return _load_gr(file_name)
    catch e
        @error "File not found"
    end
end