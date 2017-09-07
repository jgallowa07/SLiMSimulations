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
col8 = []
col9 = []
col10 = []
col11 = []
col12 = []

x_axis = []
File = open("../Output1/MyRecipe2_1_3/AveragePheno.txt","r")
first = File.readline()
numAxis,interval = (int(i) for i in first.split())
for line in File:
	n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12= (float(s) for s in line.split())
	col1.append(n1)
	col2.append(n2)
	col3.append(n3)
	col4.append(n4)
	col5.append(n5)
	col6.append(n6)
	col7.append(n7)
	col8.append(n8)
	col9.append(n9)
	col10.append(n10)
	col11.append(n11)
	col12.append(n12)
	#col3.append(n3)
	#col4.append(n4)
	#trace5.append(n5)
		
		

print("numAxis: "+str(numAxis))	

for i in range(0,numAxis + 1):
	x_axis.append(i*interval)
print(x_axis)	

trace0 = go.Scatter(
    x = x_axis,
    y = col1,
    name = 'Oceanic Average Phenotype',
    line = dict(
        color = ('rgb(25, 12, 24)'),
        width = 2)
)
trace1 = go.Scatter(
    x = x_axis,
    y = col2,
    name = 'Freshwater average phenotype',
    line = dict(
        color = ('rgb(225, 96, 167)'),
        width = 2,)
)
trace2 = go.Scatter(
    x = x_axis,
    y = col3,
    name = 'Lake 1 average phenotype',
    line = dict(
        color = ('rgb(5, 16, 117)'),
        width = 1,)
)
trace3 = go.Scatter(
    x = x_axis,
    y = col4,
    name = 'Lake 2 average phenotype',
    line = dict(
        color = ('rgb(25, 196, 67)'),
        width = 1,)
)
trace4 = go.Scatter(
    x = x_axis,
    y = col5,
    name = 'Lake 3 average phenotype',
    line = dict(
        color = ('rgb(25, 196, 67)'),
        width = 1,)
)
trace5 = go.Scatter(
    x = x_axis,
    y = col6,
    name = 'Lake 4 average phenotype',
    line = dict(
        color = ('rgb(25, 196, 67)'),
        width = 1,)
)
trace6 = go.Scatter(
    x = x_axis,
    y = col7,
    name = 'Lake 5 average phenotype',
    line = dict(
        color = ('rgb(25, 196, 67)'),
        width = 1,)
)
trace7 = go.Scatter(
    x = x_axis,
    y = col8,
    name = 'Lake 6 average phenotype',
    line = dict(
        color = ('rgb(25, 196, 67)'),
        width = 1,)
)
trace8 = go.Scatter(
    x = x_axis,
    y = col9,
    name = 'Lake 7 average phenotype',
    line = dict(
        color = ('rgb(25, 196, 67)'),
        width = 1,)
)
trace9 = go.Scatter(
    x = x_axis,
    y = col10,
    name = 'Lake 8 average phenotype',
    line = dict(
        color = ('rgb(25, 196, 67)'),
        width = 1,)
)
trace10 = go.Scatter(
    x = x_axis,
    y = col11,
    name = 'Lake 9 average phenotype',
    line = dict(
        color = ('rgb(25, 196, 67)'),
        width = 1,)
)
trace11 = go.Scatter(
    x = x_axis,
    y = col12,
    name = 'Lake 10 average phenotype',
    line = dict(
        color = ('rgb(25, 196, 67)'),
        width = 1,)
)
print(type(trace11))
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
data = [trace0, trace1, trace2, trace3,trace4,trace5,trace6,trace7,trace8,trace9,trace10,trace11 ]

# Edit the layout
layout = dict(title = 'Average Phenotype',
              xaxis = dict(title = 'Generations'),
              yaxis = dict(title = 'Phenotype'),
              )

fig = dict(data=data, layout=layout)
py.iplot(fig, filename='2_1_3 Average Phenotype')
