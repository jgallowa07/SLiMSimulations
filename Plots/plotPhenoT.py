import plotly
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go
import plotly.plotly as py
import sys




data = []
x_axis = []
lake_traces = []
File = open("../Output1/MyRecipe2_2_9/AveragePheno.txt","r")
first = File.readline()
numAxis,interval,numLakes = (int(i) for i in first.split())

for i in range(0,(numLakes+2)):
	data.append([])

for line in File:
	count = 0
	for s in line.split():
		data[count].append(float(s))
		count += 1
		
print("numAxis: "+str(numAxis))	

for i in range(0,numAxis + 1):
	x_axis.append(i*interval)
print(x_axis)	



trace0 = go.Scatter(
    x = x_axis,
    y = data[0],
    name = 'Marine Average Phenotype',
    line = dict(
        color = ('rgb(25, 12, 24)'),
        width = 2)
)
trace1 = go.Scatter(
    x = x_axis,
    y = data[1],
    name = 'Freshwater average Phenotype',
    line = dict(
        color = ('rgb(225, 96, 167)'),
        width = 2,)
)

for i in range(0,numLakes):

	numlake = i+1	

	trace = go.Scatter(
		x = x_axis,
		y = data[i+2],	
		name = 'lake '+str(numlake)+ ' average Phenotype',
		line = dict(
			color = ('rgb(25,196,67)'),
			width = 1,)
	)
	
	lake_traces.append(trace)	


data = [trace0, trace1, lake_traces[0], lake_traces[3] , lake_traces[6], lake_traces[9]]

# Edit the layout
layout = dict(title = 'Average Phenotype',
              xaxis = dict(title = 'Generations'),
              yaxis = dict(title = 'Phenotype'),
              )

fig = dict(data=data, layout=layout)
py.iplot(fig, filename='2_2_9 Average Phenotype')
