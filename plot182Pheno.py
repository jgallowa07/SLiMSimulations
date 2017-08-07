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

x_axis = []
numAxis = 40

def main():
	File = open("myOut178_Ph","r")
	for line in File:
		n1, n2= (float(s) for s in line.split())
		col1.append(n1)
		col2.append(n2)
		#col3.append(n3)
		#col4.append(n4)
		#trace5.append(n5)
		
	
	return 1

		
main();


#print("numAxis: "+str(numAxis))	

for i in range(0,numAxis + 1):
	x_axis.append(i*500)
print(x_axis)	

trace0 = go.Scatter(
    x = x_axis,
    y = col1,
    name = 'Oceanic Average Phenotype',
    line = dict(
        color = ('rgb(25, 12, 24)'),
        width = 1)
)
trace1 = go.Scatter(
    x = x_axis,
    y = col2,
    name = 'Freshwater average phenotype',
    line = dict(
        color = ('rgb(225, 96, 167)'),
        width = 1,)
)
#trace2 = go.Scatter(
#    x = x_axis,
#    y = col3,
#    name = 'High Frequency m2',
#    line = dict(
#	color = ('rgb(5, 125, 24)'),
#	width = 1) # dash options include 'dash', 'dot', and 'dashdot'
#)
#trace3 = go.Scatter(
#    x = x_axis,
#    y = col4,
#    name = 'Low Frequency m2',
#    line = dict(
#	color = ('rgb(22, 96, 167)'),
#	width = 1)
#)
#trace4 = go.Scatter(
#    x = month,
#    y = high_2000,
#    name = 'High 2000',
#    line = dict(
#        color = ('rgb(205, 12, 24)'),
#        width = 1,
#        dash = 'dot')
#)
#trace5 = go.Scatter(
#    x = month,
#    y = low_2000,
#    name = 'Low 2000',
#    line = dict(
#        color = ('rgb(22, 96, 167)'),
#        width = 1,
#        dash = 'dot')
#)
data = [trace0, trace1]

# Edit the layout
layout = dict(title = 'Average Phenotype',
              xaxis = dict(title = 'Generations'),
              yaxis = dict(title = 'Phenotype'),
              )

fig = dict(data=data, layout=layout)
py.iplot(fig, filename='MyRecipe1_7_8(Phenotypes)')
