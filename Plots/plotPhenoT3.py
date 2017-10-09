import plotly
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go
import plotly.plotly as py
import sys

col1 = []
col2 = []
col3 = []


data = []
x_axis = []
x2_axis = []
lake_traces = []
File = open("../Output1/MyRecipe4_0_3Local/AveragePheno.txt","r")
first = File.readline()
numAxis,interval,newLakesIntroduced = (int(i) for i in first.split())

#for i in range(0,(numLakes+2)):
#	data.append([])

count = 0
for line in File:
	if(count < newLakesIntroduced):
		print("made it")
		n1, n2= (float(s) for s in line.split())
		col1.append(n1)
		col2.append(n2)
	else:
		n1, n2, n3= (float(s) for s in line.split())
		col1.append(n1)
		col2.append(n2)
		col3.append(n3)
	count += 1
		
print("numAxis: "+str(numAxis))	

for i in range(1,numAxis + 1):
	x_axis.append(i*interval)
print(x_axis)	

x2_axis = x_axis[newLakesIntroduced:]

trace0 = go.Scatter(
    x = x_axis,
    y = col1,
    name = 'Marine Average Phenotype',
    line = dict(
        color = ('rgb(25, 12, 24)'),
        width = 2)
)
trace1 = go.Scatter(
    x = x_axis,
    y = col2,
    name = 'Freshwater average Phenotype',
    line = dict(
        color = ('rgb(225, 96, 167)'),
        width = 2,)
)

trace2 = go.Scatter(
	x = x2_axis,
	y = col3,	
	name = 'new Lakes average Phenotype',
	line = dict(
		color = ('rgb(25,196,67)'),
		width = 2,)
)
	
data = [trace0, trace1, trace2]

# Edit the layout
layout = dict(title = 'Average Phenotype',
              xaxis = dict(title = 'Generations'),
              yaxis = dict(title = 'Phenotype'),
              )

fig = dict(data=data, layout=layout)
py.iplot(fig, filename='4_0_3 Average Phenotype')