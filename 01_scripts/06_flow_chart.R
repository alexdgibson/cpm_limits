# flow_chart.R
# create a flow chart for the analysis and screening process including pilot analysis
# January 2025

# load in the DiagrammeR package
library(DiagrammeR)


# create the graph
grViz("digraph {
  graph[layout = dot, 
        rankdir = TB,
        overlap = true,
        fontsize = 10]
  node [shape = rectangle,
  fixedsize = true,
  width = 4.2]
  
  'Develop statistical code/methods (EuroSCORE II)'
  'Test statistical code/methods (Framingham)'
  'Complete search'
  'Screen articles'
  'Collect data'
  'Pilot analysis on neurology'
  'Potential protocol changes'
  'Complete final analysis'
  
  'Develop statistical code/methods (EuroSCORE II)' -> 'Test statistical code/methods (Framingham)' ->
  'Complete search' ->'Screen articles' ->'Collect data' ->'Pilot analysis on neurology' ->'Potential protocol changes' ->
  'Complete final analysis'

      }")

