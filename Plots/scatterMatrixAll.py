import plotly.graph_objs as go
import plotly.plotly as py

import sys

col1 = []
col2 = []

File = open("../Output1/MyRecipe3_0_1/MatrixAllelesAll.txt","r")
first = File.readline()
for line in File:
	n1, n2 = (float(s) for s in line.split())
	col1.append(n1)
	col2.append(n2)

trace1 = go.Scatter(
    x = col2,
    y = col1,
    mode='markers',
    marker=dict(
        size='3',
        color = '#FFBAD2'
    )
)


data = [trace1]

layout = dict(title = 'Sample Haplotypes',
              xaxis = dict(title = 'Loci Positions'),
              yaxis = dict(title = 'Genomes'),
              )

fig = dict(data = data, layout = layout)
	
py.iplot(fig, filename='3_0_1 Matrix Alelles All Test')
