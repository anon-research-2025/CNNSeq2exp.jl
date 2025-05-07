using Printf
const markers = [:+, :diamond, :star4, :rtriangle, :rect,
    :circle, :pentagon, :cross, :star5]


function FigGridScatters(
        predicted::Vector{Vector{T}}, labeled::Vector{Vector{T}}, 
        titles::Vector{String}, r_values::Vector{T}, r2_values::Vector{T};
        fig_size = (1800, 350),
        markersize=8, strokewidth=1
        ) where T <: Real
    len_predicted = length(predicted)
    @assert len_predicted == length(labeled) "Lengths of predicted and labeled must be equal"
    fig = Figure(size = fig_size, )
        
        for i in 1:len_predicted
            ax = Axis(fig[1, i], aspect = AxisAspect(1), backgroundcolor=:transparent)
            CairoMakie.scatter!(ax, predicted[i], labeled[i];  
                markersize = markersize, marker = :circle, strokewidth = strokewidth, 
                strokecolor = :black, color = :orange,)
            limits!(ax, -1.1, 1.1, -1.45, 1.1) # assume tanh output
            ax.xticks = [-1, 0, 1]
            ax.yticks = [-1, 0, 1]

            if i != 1
                ax.yticksvisible = false
                ax.yticklabelsvisible = false
            end

            ax.xticklabelsize = 32
            ax.yticklabelsize = 32
            ax.title = titles[i]
            ax.titlealign = :center
            ax.titlefont = :bold
            ax.titlesize = 37
            ax.titlecolor = :black
            ax.titlegap = 10
            ax.xlabelsize = 28
            ax.ylabelsize = 28
            
            # Labels
            ax.xlabel = "Predicted values"
            if i == 1
                ax.ylabel = "Labeled values"  # Only on the leftmost plot
            end
        
            # Corrected text placement with `rich()`
            text!(ax, [-0.85], [-1.2], 
                text = rich(rich("r = "), @sprintf("%.2f", r_values[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
            text!(ax, [0.15], [-1.2], 
                text = rich(rich("R² = "), @sprintf("%.2f", r2_values[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
        end
    colgap!(fig.layout, 0.25)
    fig
end


# use case
# fig = FigGridScatters(predicted, labeled, titles, r_values, r2_values;
#     fig_size = (2400, 450)
#     )

# save("pub_figs/scatter_plots.png", fig, px_per_unit=1)


function FigGridScatters2(
        predicted::Vector{Vector{T}}, 
        predicted2::Vector{Vector{T}}, 
        labeled::Vector{Vector{T}}, 
        labeled2::Vector{Vector{T}}, 
        titles::Vector{String}, 
        r_values::Vector{T}, 
        r_values2::Vector{T}, 
        r2_values::Vector{T},
        r2_values2::Vector{T};
        fig_size = (1800, 350),
        markersize=8, strokewidth=1
        ) where T <: Real
    len_predicted = length(predicted)
    @assert len_predicted == length(labeled) "Lengths of predicted and labeled must be equal"
    fig = Figure(size = (fig_size), )
        
        for i in 1:len_predicted
            ax = Axis(fig[1, i], aspect = AxisAspect(1), backgroundcolor=:transparent)
            CairoMakie.scatter!(ax, predicted[i], labeled[i];  
                markersize = markersize, marker = :circle, strokewidth = strokewidth, 
                strokecolor = :black, color = :orange,)
            limits!(ax, -1.1, 1.1, -1.45, 1.1) # assume tanh output
            ax.xticks = [-1, 0, 1]
            ax.yticks = [-1, 0, 1]

            if i != 1
                ax.yticksvisible = false
                ax.yticklabelsvisible = false
            end

            ax.xticklabelsize = 32
            ax.yticklabelsize = 32
            ax.title = titles[i]
            ax.titlealign = :center
            ax.titlefont = :bold
            ax.titlesize = 28
            ax.titlecolor = :black
            ax.titlegap = 10
            ax.xlabelsize = 28
            ax.ylabelsize = 28
            ax.xticksvisible = false
            ax.xticklabelsvisible = false
            ax.xlabelvisible = false
            
            # Labels
            ax.xlabel = "Predicted values"
            if i == 1
                ax.ylabel = "Labeled values"  # Only on the leftmost plot
            end
        
            # Corrected text placement with `rich()`
            text!(ax, [-0.85], [-1.2], 
                text = rich(rich("r = "), @sprintf("%.2f", r_values[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
            text!(ax, [0.15], [-1.2], 
                text = rich(rich("R² = "), @sprintf("%.2f", r2_values[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
        end

        for i in 1:len_predicted
            ax = Axis(fig[2, i], aspect = AxisAspect(1), backgroundcolor=:transparent)
            CairoMakie.scatter!(ax, predicted2[i], labeled2[i];  
                markersize = markersize, marker = :circle, strokewidth = strokewidth, 
                strokecolor = :black, color = :orange,)
            limits!(ax, -1.1, 1.1, -1.45, 1.1) # assume tanh output
            ax.xticks = [-1, 0, 1]
            ax.yticks = [-1, 0, 1]

            if i != 1
                ax.yticksvisible = false
                ax.yticklabelsvisible = false
            end

            ax.xticklabelsize = 32
            ax.yticklabelsize = 32
            # ax.title = titles[i]
            ax.titlealign = :center
            ax.titlefont = :bold
            ax.titlesize = 28
            ax.titlecolor = :black
            ax.titlegap = 10
            ax.xlabelsize = 28
            ax.ylabelsize = 28
            
            # Labels
            ax.xlabel = "Predicted values"
            if i == 1
                ax.ylabel = "Labeled values"  # Only on the leftmost plot
            end
        
            # Corrected text placement with `rich()`
            text!(ax, [-0.85], [-1.2], 
                text = rich(rich("r = "), @sprintf("%.2f", r_values2[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
            text!(ax, [0.15], [-1.2], 
                text = rich(rich("R² = "), @sprintf("%.2f", r2_values2[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
        end

        colgap!(fig.layout, 15)
        rowgap!(fig.layout, 50)

        o = 25
        Label(fig.layout[1, :, TopLeft()], "A",
            fontsize = 41,
            font = :bold,
            padding = (0, o, o, 0),
            halign = :right)
        
        Label(fig.layout[2, :, TopLeft()], "B",
            fontsize = 41,
            font = :bold,
            padding = (0, o, o, 0),
            halign = :right)
    fig
end


function FigGridScatters3(
        predicted::Vector{Vector{T}}, 
        predicted2::Vector{Vector{T}}, 
        predicted3::Vector{Vector{T}},
        labeled::Vector{Vector{T}}, 
        labeled2::Vector{Vector{T}}, 
        labeled3::Vector{Vector{T}},
        titles::Vector{String}, 
        r_values::Vector{T}, 
        r_values2::Vector{T}, 
        r_values3::Vector{T},
        r2_values::Vector{T},
        r2_values2::Vector{T},
        r2_values3::Vector{T};
        fig_size = (1800, 350),
        markersize=8, strokewidth=1
        ) where T <: Real
    len_predicted = length(predicted)
    @assert len_predicted == length(labeled) "Lengths of predicted and labeled must be equal"
    fig = Figure(size = (fig_size), )
        
        for i in 1:len_predicted
            ax = Axis(fig[1, i], aspect = AxisAspect(1), backgroundcolor=:transparent)
            CairoMakie.scatter!(ax, predicted[i], labeled[i];  
                markersize = markersize, marker = :circle, strokewidth = strokewidth, 
                strokecolor = :black, color = :orange,)
            limits!(ax, -1.1, 1.1, -1.45, 1.1) # assume tanh output
            ax.xticks = [-1, 0, 1]
            ax.yticks = [-1, 0, 1]

            if i != 1
                ax.yticksvisible = false
                ax.yticklabelsvisible = false
            end

            ax.xticklabelsize = 32
            ax.yticklabelsize = 32
            ax.title = titles[i]
            ax.titlealign = :center
            ax.titlefont = :bold
            ax.titlesize = 28
            ax.titlecolor = :black
            ax.titlegap = 10
            ax.xlabelsize = 28
            ax.ylabelsize = 28
            ax.xticksvisible = false
            ax.xticklabelsvisible = false
            ax.xlabelvisible = false
            
            # Labels
            ax.xlabel = "Predicted values"
            if i == 1
                ax.ylabel = "Labeled values"  # Only on the leftmost plot
            end
        
            # Corrected text placement with `rich()`
            text!(ax, [-0.85], [-1.2], 
                text = rich(rich("r = "), @sprintf("%.2f", r_values[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
            text!(ax, [0.15], [-1.2], 
                text = rich(rich("R² = "), @sprintf("%.2f", r2_values[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
        end

        for i in 1:len_predicted
            ax = Axis(fig[2, i], aspect = AxisAspect(1), backgroundcolor=:transparent)
            CairoMakie.scatter!(ax, predicted2[i], labeled2[i];  
                markersize = markersize, marker = :circle, strokewidth = strokewidth, 
                strokecolor = :black, color = :orange,)
            limits!(ax, -1.1, 1.1, -1.45, 1.1) # assume tanh output
            ax.xticks = [-1, 0, 1]
            ax.yticks = [-1, 0, 1]

            if i != 1
                ax.yticksvisible = false
                ax.yticklabelsvisible = false
            end

            ax.xticklabelsize = 32
            ax.yticklabelsize = 32
            # ax.title = titles[i]
            ax.titlealign = :center
            ax.titlefont = :bold
            ax.titlesize = 28
            ax.titlecolor = :black
            ax.titlegap = 10
            ax.xlabelsize = 28
            ax.ylabelsize = 28
            ax.xticksvisible = false
            ax.xticklabelsvisible = false
            ax.xlabelvisible = false
            
            
            # Labels
            ax.xlabel = "Predicted values"
            if i == 1
                ax.ylabel = "Labeled values"  # Only on the leftmost plot
            end
        
            # Corrected text placement with `rich()`
            text!(ax, [-0.85], [-1.2], 
                text = rich(rich("r = "), @sprintf("%.2f", r_values2[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
            text!(ax, [0.15], [-1.2], 
                text = rich(rich("R² = "), @sprintf("%.2f", r2_values2[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
        end

        for i in 1:len_predicted
            ax = Axis(fig[3, i], aspect = AxisAspect(1), backgroundcolor=:transparent)
            CairoMakie.scatter!(ax, predicted3[i], labeled3[i];  
                markersize = markersize, marker = :circle, strokewidth = strokewidth, 
                strokecolor = :black, color = :orange,)
            limits!(ax, -1.1, 1.1, -1.45, 1.1) # assume tanh output
            ax.xticks = [-1, 0, 1]
            ax.yticks = [-1, 0, 1]

            if i != 1
                ax.yticksvisible = false
                ax.yticklabelsvisible = false
            end

            ax.xticklabelsize = 32
            ax.yticklabelsize = 32
            # ax.title = titles[i]
            ax.titlealign = :center
            ax.titlefont = :bold
            ax.titlesize = 28
            ax.titlecolor = :black
            ax.titlegap = 10
            ax.xlabelsize = 28
            ax.ylabelsize = 28
            
            # Labels
            ax.xlabel = "Predicted values"
            if i == 1
                ax.ylabel = "Labeled values"  # Only on the leftmost plot
            end
        
            # Corrected text placement with `rich()`
            text!(ax, [-0.85], [-1.2], 
                text = rich(rich("r = "), @sprintf("%.2f", r_values3[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
            text!(ax, [0.15], [-1.2], 
                text = rich(rich("R² = "), @sprintf("%.2f", r2_values3[i]), font=:bold), 
                fontsize = 25, align = (:left, :center), color = :black
            )
        end

        colgap!(fig.layout, 15)
        rowgap!(fig.layout, 35)

        o = 20
        fs = 36
        Label(fig.layout[1, :, Top()], "CNN",
            fontsize = fs,
            font = :bold,
            padding = (0, o, o, 0),
            halign = :left)
        
        Label(fig.layout[2, :, Top()], "Linear regression: additive",
            fontsize = fs,
            font = :bold,
            padding = (0, o, o, 0),
            halign = :left)

        Label(fig.layout[3, :, Top()], "Linear regression: paired",
            fontsize = fs,
            font = :bold,
            padding = (0, o, o, 0),
            halign = :left)
    fig
end