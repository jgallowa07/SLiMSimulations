import plotly
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go
import plotly.plotly as py
import sys

dataFromFile = []
dataFromFile2 = []
traceNames = []
traceColors = []
traceNames2 = []
traceColors2 = []
xAxis = []
xAxis2 = []
Traces = []

File = open("../Output1/MyRecipe8_0_6/AverageFst.txt","r")
File2 = open("../Output1/MyRecipe8_0_6/AverageFst2.txt","r")
first = File.readline()
first2 = File2.readline()
numAxis,interval = (int(i) for i in first.split())
numAxis2,interval2,introduce = (int(i) for i in first2.split())

for i in range(0,9):
	col = []
	dataFromFile.append(col)

for i in range(0,9):
	col = []
	dataFromFile2.append(col)

for line in File:
	lineData = line.split()
	for i in range(0,9):
		dataFromFile[i].append(lineData[i])

for line in File2:
	lineData = line.split()
	for i in range(0,9):
		dataFromFile2[i].append(lineData[i])

for i in range(0,numAxis + 1):
	xAxis.append(i*interval)

for i in range(0,numAxis2 + 1):
	xAxis2.append(i*interval2+introduce)	


trace1name = 'Marine,Lakes'
trace1color = 'rgb(178,34,34)'
traceNames.append(trace1name)
traceColors.append(trace1color)

trace2name = 'Marine,Lakes:Neutral'
trace2color = 'rgb(233,150,122)'
traceNames.append(trace2name)
traceColors.append(trace2color)

trace3name = 'Marine,Lakes:Effect'
trace3color = 'rgb(255,140,0)'
traceNames.append(trace3name)
traceColors.append(trace3color)

trace4name = 'Lakes'
trace4color = 'rgb(189,183,107)'
traceNames.append(trace4name)
traceColors.append(trace4color)

trace5name = 'Lakes:Neutral'
trace5color = 'rgb(128,128,0)'
traceNames.append(trace5name)
traceColors.append(trace5color)

trace6name = 'Lakes:Effect'
trace6color = 'rgb(34,139,34)'
traceNames.append(trace6name)
traceColors.append(trace6color)

trace7name = 'Marine'
trace7color = 'rgb(47,27,27)'
traceNames.append(trace7name)
traceColors.append(trace7color)

trace7name = 'Marine:Neutral'
trace7color = 'rgb(147,27,27)'
traceNames.append(trace7name)
traceColors.append(trace7color)

trace7name = 'Marine:Effect'
trace7color = 'rgb(47,127,27)'
traceNames.append(trace7name)
traceColors.append(trace7color)

#--------------------------

trace8name = 'Marine,Introduced'
trace8color = 'rgb(0,128,128)'
traceNames2.append(trace8name)
traceColors2.append(trace8color)

trace9name = 'Marine,Intro:Neutral'
trace9color = 'rgb(64,224,208)'
traceNames2.append(trace9name)
traceColors2.append(trace9color)

trace10name = 'Marine,Intro:Effect'
trace10color = 'rgb(127,255,212)'
traceNames2.append(trace10name)
traceColors2.append(trace10color)

trace11name = 'Introduced'
trace11color = 'rgb(139,0,139)'
traceNames2.append(trace11name)
traceColors2.append(trace11color)

trace12name = 'Introduced:Neutral'
trace12color = 'rgb(216,160,221)'
traceNames2.append(trace12name)
traceColors2.append(trace12color)

trace13name = 'Introduced:Effect'
trace13color = 'rgb(199,21,133)'
traceNames2.append(trace13name)
traceColors2.append(trace13color)

trace14name = 'Lakes,Introduced'
trace14color = 'rgb(160,82,45)'
traceNames2.append(trace14name)
traceColors2.append(trace14color)

trace15name = 'Lakes,Intro:Neutral'
trace15color = 'rgb(210,105,30)'
traceNames2.append(trace15name)
traceColors2.append(trace15color)

trace16name = 'Lakes,Intro:Effect'
trace16color = 'rgb(188,143,143)'
traceNames2.append(trace16name)
traceColors2.append(trace16color)


for i in range(0,9):

	trace = go.Scatter(
	    x = xAxis,
	    y = dataFromFile[i],
	    name = traceNames[i],
	    line = dict(
		color = (traceColors[i]),
		width = 1)
	)

	Traces.append(trace)
	
for i in range(0,9):

	trace = go.Scatter(
	    x = xAxis2,
	    y = dataFromFile2[i],
	    name = traceNames2[i],
	    line = dict(
		color = (traceColors2[i]),
		width = 1)
	)

	Traces.append(trace)

data = Traces

# Edit the layout
layout = dict(title = 'Average Fst',
              xaxis = dict(title = 'Generations'),
              yaxis = dict(title = 'Fst'),
              )

fig = dict(data=data, layout=layout)
#py.iplot(fig, filename='3_0_1 Average Fst')
plotly.offline.plot(fig, filename='8_0_6Fst.html')
