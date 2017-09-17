import sys
import numpy as np
import plotly.plotly as py
import plotly.graph_objs as go

#programmers = ['Alex','Nicole','Sara','Etienne','Chelsea','Jody','Marianne']

#base = datetime.datetime.today()
#date_list = [base - datetime.timedelta(days=x) for x in range(0, 180)]
x = []
y = []
z = []

File = open("../Output1/MyRecipe3_0_0Local/MatrixAlleles.txt","r")
first = File.readline()
numGenomes,numLoci = (int(i) for i in first.split())
	
for i in range(0,numLoci):
	x.append(i)

for i in range(0,numGenomes):
	y.append(i)
	
for line in File:
	newRow = []
	for pos in line.split():
		newRow.append(pos)
	z.append(list(newRow))
	

#for prgmr in programmers:
#    new_row = []
#    for date in date_list:
#        new_row.append( np.random.poisson() )
#    z.append(list(new_row))

data = [
    go.Heatmap(
        z=z,
        x=x,
        y=y,
        colorscale='Viridis',
    )
]

layout = go.Layout(
    title='Alleles',
    xaxis = dict(ticks='', nticks=36),
    yaxis = dict(ticks='' )
)

fig = go.Figure(data=data, layout=layout)
py.iplot(fig, filename='Test HeatMap')
