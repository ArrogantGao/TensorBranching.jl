using Graphs, TensorBranching
using CairoMakie
using GenericTensorNetworks
using TensorBranching.OMEinsum, TensorBranching.GenericTensorNetworks, TensorBranching.AbstractTrees
using Statistics, LsqFit
using LaTeXStrings

function three_rr_line(n_in::Int, N_sub::Int)
    g = SimpleGraph(n_in * N_sub)
    for i in 1:N_sub
        g_i = random_regular_graph(n_in, 3)
        for edge in edges(g_i)
            add_edge!(g, (i-1)*n_in + src(edge), (i-1)*n_in + dst(edge))
        end
        (i > 1) && add_edge!(g, (i - 1) * n_in, (i - 1) * n_in + 1)
    end
    return g
end

begin
    n_in = 10
    N_subs = [1:15...]
    branch_num = []
    cc_num = []
    for N_sub in N_subs
        count_i = []
        cc_i = []
        for i in 1:20
            g = three_rr_line(n_in, N_sub)
            count_mis2 = counting_mis2(g)
            push!(count_i, count_mis2.count)

            prob = GenericTensorNetwork(IndependentSet(g); optimizer=TreeSA(; ntrials=1, niters=10))
            code = prob.code.eins

            cc = contraction_complexity(code, uniformsize(code, 2))
            push!(cc_i, cc.tc)

            @show N_sub, i, count_mis2.count, cc.tc
        end
        push!(branch_num, mean(count_i))
        push!(cc_num, mean(cc_i))
    end
end

begin    
    @. model_branch(x, p) = p[1]^x
    fit_branch = curve_fit(model_branch, n_in .* N_subs, branch_num, [1.0])
    @show fit_branch.param

    @. model_tn(x, p) = p[1] * x^p[2]
    fit_tn = curve_fit(model_tn, n_in .* N_subs, 2.0 .^ cc_num, [1.0, 1.0])
    @show fit_tn.param
end


begin
    f = Figure()
    ax1 = Axis(f[1, 1], yticklabelcolor = :blue, yscale = log10, xlabel = "Number of sub-networks", ylabel = "Number of branches")
    ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right, yscale = log2, ylabel = "Contraction complexity")
    hidespines!(ax2)
    hidexdecorations!(ax2)

    scatter!(ax1, N_subs, branch_num, markersize = 10, label = "branching", color = :blue)
    scatter!(ax2, N_subs, 2.0 .^ cc_num, markersize = 10, label = "tn", color = :red)

    lines!(ax1, N_subs, model_branch(n_in .* N_subs, fit_branch.param), label = "fit", color = :blue)
    lines!(ax2, N_subs, model_tn(n_in .* N_subs, fit_tn.param), label = "fit", color = :red)

    text!(ax1, 8, 100, text = L"$O({%$(round(fit_branch.param[1], digits = 4))}^n)$", color = :blue, fontsize = 20)
    text!(ax2, 5, 2^11, text = L"$O({n^{%$(round(fit_tn.param[2], digits = 2))}})$", color = :red, fontsize = 20)
end
f

save("numbranch.pdf", f)