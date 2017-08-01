import plotly
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go
import plotly.plotly as py
import sys
trace1 = []
trace2 = []
trace3 = []
trace4 = []
trace5 = []

x_axis = []
numAxis = 0

def main():
	File = open("myOut144.txt","r")
	for i in range(0,15):
		File.readline()
	numAxis = int(File.readline())
	print("numAxis: "+str(numAxis))	
	for line in File:
		n1, n2, n3, n4 = (float(s) for s in line.split())
		trace1.append(n1)
		trace2.append(n2)
		trace3.append(n3)
		trace4.append(n4)
		#trace5.append(n5)
		
	
	return numAxis

		
numAxis = main();


print("numAxis: "+str(numAxis))	

for i in range(0,numAxis + 1):
	x_axis.append(i*50)
print(x_axis)	

Trace1 = go.Scatter(x = x_axis,y = trace1)
Trace2 = go.Scatter(x = x_axis,y = trace2)
Trace3 = go.Scatter(x = x_axis,y = trace3)
Trace4 = go.Scatter(x = x_axis,y = trace4)
#Trace5 = go.Scatter(x = x_axis,y = trace5)

py.iplot([Trace1,Trace2,Trace3,Trace4],filename='1_4_4')


#plotly.offline.plot({
#    "data": [Scatter(x=[1, 2, 3, 4], y=[4, 3, 2, 1.1])],
#    "layout": Layout(title="hello world")
#})
