"""
Plot utility functions

Author: Jonathan Richmond
C: 6/3/25
"""

using MBD, GLMakie, LaTeXStrings

export set2DPlotParameters, set3DPlotParameters

"""
    primary_resolution()

Return screen dimensions
"""
function primary_resolution()
    monitor = GLMakie.GLFW.GetPrimaryMonitor()
    videoMode = GLMakie.MonitorProperties(monitor).videomode
    return (videoMode.width, videoMode.height)
end

"""
    set2DPlotParameters(plotTitle, plotXLabel, plotYLabel)

Return 2-D plot figure and axes

# Arguments
- `plotTitle::LaTeXString`: Plot title
- `plotXLabel::LaTeXString`: Plot x-axis label
- `plotYLabel::LaTeXString`: Plot y-axis label
"""
function set2DPlotParameters(plotTitle::LaTeXStrings.LaTeXString, plotXLabel::LaTeXStrings.LaTeXString, plotYLabel::LaTeXStrings.LaTeXString)
    GLMakie.set_theme!(GLMakie.theme_black())
    figure = GLMakie.Figure(size = primary_resolution(), font = "CMU Serif", fontsize = 30)
    axis = GLMakie.Axis(figure[1,1], aspect = GLMakie.DataAspect(), title = plotTitle, xlabel = plotXLabel, ylabel = plotYLabel)

    return (figure, axis)
end

"""
    set3DPlotParameters(plotTitle, plotXLabel, plotYLabel, plotZLabel)

Return 3-D plot figure and axes

# Arguments
- `plotTitle::LaTeXString`: Plot title
- `plotXLabel::LaTeXString`: Plot x-axis label
- `plotYLabel::LaTeXString`: Plot y-axis label
- `plotZLabel::LaTeXString`: Plot z-axis label
"""
function set3DPlotParameters(plotTitle::LaTeXStrings.LaTeXString, plotXLabel::LaTeXStrings.LaTeXString, plotYLabel::LaTeXStrings.LaTeXString, plotZLabel::LaTeXStrings.LaTeXString)
    GLMakie.set_theme!(GLMakie.theme_black())
    figure = GLMakie.Figure(size = primary_resolution(), font = "CMU Serif", fontsize = 30)
    axis = GLMakie.Axis3(figure[1,1], aspect = GLMakie.DataAspect(), title = plotTitle, xlabel = plotXLabel, ylabel = plotYLabel, zlabel = plotZLabel)

    return (figure, axis)
end
