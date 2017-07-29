import plotly.plotly as py
import plotly.graph_objs as go

import numpy as np

col1 = []
col2 = []
col3 = []
col4 = []
col5 = []

x_axis = []
numAxis = 0

def main():
	File = open("myOut163.txt","r")
	for i in range(0,15):
		File.readline()
	numAxis = int(File.readline())
	print("numAxis: "+str(numAxis))	
	for line in File:
		#n1 = (float(s) for s in line.split())
		n1 = int(line)
		col1.append(n1)
		#col2.append(n2)
		#col3.append(n3)
		#col4.append(n4)
		#trace5.append(n5)
		
	
	return numAxis

		
numAxis = main();
print(col1)
#x = np.random.randn(500)
data = [go.Histogram(x=col1)]

py.iplot(data, filename='basic histogram')
