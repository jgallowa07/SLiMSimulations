import plotly.graph_objs as go
import plotly.plotly as py

import sys

col1 = []
col2 = []

File = open("../Output1/MyRecipe2_2_9/OceanFreshwaterFst.txt","r")
first = File.readline()
for line in File:
	n1, n2 = (float(s) for s in line.split())
	col1.append(n1)
	col2.append(n2)

trace1 = go.Scatter(
    x = col2,
    y = col1,
    mode='markers',
    marker=dict(
        size='8',
        color = '#FFBAD2'
    )
)


data = [trace1]

layout = dict(title = 'Effect and Frequecy by Position',
              xaxis = dict(title = 'Position'),
              yaxis = dict(title = 'Fst'),
              )

fig = dict(data = data, layout = layout)
	
py.iplot(fig, filename='2_2_9 Ocean-Freshwater Fst')
