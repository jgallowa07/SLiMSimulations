import plotly.plotly as py
import plotly
import plotly.graph_objs as go

import numpy as np

col1 = []


File = open("../Output1/MyRecipe5_0_1/PhenoDist.txt","r")
#First = File.readline()	
	
for line in File:
	col1.append(float(line))
		

trace1 = go.Histogram(
    x=col1,
    histnorm='count',
    marker=dict(
        color='#6974f3',
    )
)
data = [trace1]

layout = go.Layout(
    title='Phenotype Distribution',
    xaxis=dict(
        title='Phenotype'
    ),
    yaxis=dict(
        title='# of Individuals'
    ),
    bargap=1.5,
    bargroupgap=3.2
)
fig = go.Figure(data=data, layout=layout)
#py.iplot(fig, filename='5_0_1 Phenotype Distribution')


#x = np.random.randn(500)
plotly.offline.plot(fig, filename='5_0_3 Average Phenotype')
#data = [go.Histogram(x=col1)]
#py.iplot(data, fisizeame='basic histogram')
