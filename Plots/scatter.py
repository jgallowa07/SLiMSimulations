import plotly.graph_objs as go
import plotly.plotly as py

import sys

col1 = []
col2 = []
col3 = []

File = open("../Output1/MyRecipe1_9_9/WholeMutations.txt","r")
first = File.readline()
for line in File:
	n1, n2, n3= (float(s) for s in line.split())
	col1.append(n1)
	col2.append(n2)
	col3.append(n3)

trace1 = go.Scatter(
    x = col2,
    y = col1,
    mode='markers',
    name = 'Effect',
    marker=dict(
        size='12',
        color = col3, #set color equal to a variable
        colorscale='Viridis',
        showscale=True,
    )
)


data = [trace1]

layout = dict(title = 'Effect and Frequecy by Position',
              xaxis = dict(title = 'Position'),
              yaxis = dict(title = 'Frequency'),
              )

fig = dict(data = data, layout = layout)
	
py.iplot(fig, filename='199 Whole Mutations')
