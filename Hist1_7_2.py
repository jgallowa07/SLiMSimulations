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
	File = open("myOut172.txt","r")
	for line in File:
		n1,n2,n3,n4 = (float(s) for s in line.split())
		col1.append(n1)
		col2.append(n2)
		col3.append(n3)
		col4.append(n4)
		#trace5.append(n5)
		
	
	return 1

main();
#print(len(col1))
#print(len(col2))
#print(len(col3))
#print(len(col4))
trace1 = go.Histogram(
    x=col1,
    histnorm='count',
    name='Number of "intermediate-frequency" m1 Mutations per Individual',
    #xbins=dict(
    #    start=-4.0,
    #    end=3.0,
    #    size=0.5
    #),
    marker=dict(
        color='#6974f3',
    ),
    opacity=0.75
)
trace2 = go.Histogram(
    x=col2,
    name='Number of "High/Low-frequency" m1 Mutations per Individual',
    #xbins=dict(
    #    start=-3.0,
    #    end=4,
    #    size=0.5
    #),
    marker=dict(
        color='#180913'
    ),
    opacity=0.75
)
trace3 = go.Histogram(
    x=col3,
    name='Number of "intermediate-frequency" m2 Mutations per Individual',
    #xbins=dict(
    #    start=-3.0,
    #    end=4,
    #    size=0.5
    #),
    marker=dict(
        color='#0cc7d0'
    ),
    opacity=0.75
)
trace4 = go.Histogram(
    x=col4,
    name='Number of "High/Low-frequency" m2 Mutations per Individual',
    #xbins=dict(
    #    start=-3.0,
    #    end=4,
    #    size=0.5
    #),
    marker=dict(
        color='#d11291'
    ),
    opacity=0.75
)
data = [trace1,trace2,trace3,trace4]

layout = go.Layout(
    title='One-big lots small: indiv-based metrics of realized genetic architecture, Whole Pop',
    xaxis=dict(
        title='# Mutations'
    ),
    yaxis=dict(
        title='# of Individuals'
    ),
    bargap=1.5,
    bargroupgap=3.2
)
fig = go.Figure(data=data, layout=layout)
py.iplot(fig, filename='MyRecipe1_7_2(2)')


#x = np.random.randn(500)
#data = [go.Histogram(x=col1)]
#py.iplot(data, fisizeame='basic histogram')
