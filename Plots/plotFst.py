import plotly
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go
import plotly.plotly as py
import sys
col1 = []
col2 = []
col3 = []
col4 = []
col5 = []
col6 = []
col7 = []

x_axis = []
File = open("../Output1/MyRecipe2_2_9/AverageFst.txt","r")
first = File.readline()
numAxis,interval = (int(i) for i in first.split())
for line in File:
	n1,n2,n3,n4,n5,n6,n7 = (float(s) for s in line.split())
	#n1 = float(line)
	col1.append(n1)
	col2.append(n2)
	col3.append(n3)
	col4.append(n4)
	col5.append(n5)
	col6.append(n6)
	col7.append(n7)

print("numAxis: "+str(numAxis))	

for i in range(0,numAxis + 1):
	x_axis.append(i*interval)
print(x_axis)	

trace0 = go.Scatter(
    x = x_axis,
    y = col1,
    name = 'Average Fst Marine/Freshwater',
    line = dict(
        color = ('rgb(25, 12, 24)'),
        width = 1)
)
trace1 = go.Scatter(
    x = x_axis,
    y = col2,
    name = 'Average Fst between Lakes',
    line = dict(
        color = ('rgb(225, 96, 167)'),
        width = 1,)
)
trace2 = go.Scatter(
    x = x_axis,
    y = col3,
    name = 'Average Fst between Oceanic pop',
    line = dict(
	color = ('rgb(5, 125, 24)'),
	width = 1) # dash options include 'dash', 'dot', and 'dashdot'
)
trace3 = go.Scatter(
    x = x_axis,
    y = col4,
    name = 'Average Fst Marine/Freshwater (Neutral Mutations only)',
    line = dict(
	color = ('rgb(22, 96, 167)'),
	width = 1,
	dash = 'dot')
)
trace4 = go.Scatter(
    x = x_axis,
    y = col5,
    name = 'Average Fst Marine/Freshwater (Effect Mutations only)',
    line = dict(
        color = ('rgb(205, 12, 24)'),
        width = 1,
        dash = 'dot')
)
trace5 = go.Scatter(
    x = x_axis,
    y = col6,
    name = 'Average Fst Between Lakes (Neutral Mutations only)',
    line = dict(
        color = ('rgb(5, 112, 4)'),
        width = 1,
        dash = 'dot')
)
trace6 = go.Scatter(
    x = x_axis,
    y = col7,
    name = 'Average Fst Between Lakes (Effect Mutations only)',
    line = dict(
        color = ('rgb(25, 172, 124)'),
        width = 1,
        dash = 'dot')
)
data = [trace0,trace1,trace2,trace3,trace4,trace5,trace6]

# Edit the layout
layout = dict(title = 'Average Fst',
              xaxis = dict(title = 'Generations'),
              yaxis = dict(title = 'Fst'),
              )

fig = dict(data=data, layout=layout)
py.iplot(fig, filename='2_2_9 Average Fst')
