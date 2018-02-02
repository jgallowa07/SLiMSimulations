import plotly
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go
import plotly.plotly as py
import sys

RecipeNum = sys.argv[1]
data = []

def plotThis(Fi,pop1,pop2,rgb):

	File = open(Fi,"r")
	first = File.readline()

	Fst = list(map(float,File.readline().split()))
	Pos = list(map(int,File.readline().split()))

	Fst_Sorted = [fst for _,fst in sorted(zip(Pos,Fst))]
	Pos_Sorted = sorted(Pos)

	trace = go.Scatter(
	    x = Pos_Sorted,
	    y = Fst_Sorted,
	    name = 'Fst('+pop1+','+pop2+')',
	    line = dict(
		color = (rgb),
		width = 1
			)
	)

	data.append(trace)

	# Edit the layout

plotThis("../Output1/MyRecipe"+RecipeNum+"/FreshwaterFreshwater2Fst.txt","Freshwater","IntroducedFreshwater",'rgb(60,179,113)')
plotThis("../Output1/MyRecipe"+RecipeNum+"/OceanFreshwaterFst.txt","Ocean","Freshwater",'rgb(127,255,212)')
plotThis("../Output1/MyRecipe"+RecipeNum+"/OceanFreshwater2Fst.txt","Ocean","IntroducedFreshwater",'rgb(138,43,226)')

layout = dict(
	title = 'Fst per Mutation',
	xaxis = dict(title = 'Chromosome Position'),
	yaxis = dict(title = 'Fst'),
	      )

fig = dict(data=data, layout=layout)
plotly.offline.plot(fig, filename=RecipeNum + 'FstChromosome.html')


