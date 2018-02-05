import plotly
from plotly.graph_objs import Scatter, Layout
import plotly.graph_objs as go
import plotly.plotly as py
import sys

recipeNumber = sys.argv[1]

File = open("../Output1/MyRecipe"+ recipeNumber +"/EffectMutationsThrough.txt","r")
first = File.readline()
numAxis,interval,introduce = (int(i) for i in first.split())

xAxis = []
for i in range(0,numAxis + 1):
	xAxis.append(i*interval)

intervalCount = 0

MutationFrequenciesp1 = {} 			#Key = ID : Value = [freqs]
MutationFrequenciesp2 = {} 			#Key = ID : Value = [freqs]
MutationFrequenciesp3 = {} 			#Key = ID : Value = [freqs]
MutationStartTimes = {}				#Key = ID : Value = [freqs]

lostMutationIDs = []

for line in File:
	
	ID = line.split()	
	if(line[0] == "#"):
		break
	
	MarineFreqs = File.readline().split()
	FreshwaterFreqs = File.readline().split()

	seenMutsThisGeneration = []
	for i,j in enumerate(ID):
		j = int(j)
		seenMutsThisGeneration.append(j)
		if (j not in MutationFrequenciesp1):      			#we haven't seen this mutation.
			MutationFrequenciesp1[j] = []
			MutationFrequenciesp2[j] = []
			MutationFrequenciesp3[j] = []
			MutationStartTimes[j] = [intervalCount,numAxis]
		MutationFrequenciesp1[j].append(float(MarineFreqs[i]))
		MutationFrequenciesp2[j].append(float(FreshwaterFreqs[i]))
		MutationFrequenciesp3[j].append(None)
		
	for i in MutationFrequenciesp1:
		if ((i not in seenMutsThisGeneration) & (MutationStartTimes[i][1] == numAxis)):
			MutationStartTimes[i][1] = intervalCount
	
	intervalCount += 1	

	
for line in File:
	ID = line.split()
	MarineFreqs = File.readline().split()
	FreshwaterFreqs = File.readline().split()
	IntroducedFreqs = File.readline().split()
	
	seenMutsThisGeneration = []
	for i,j in enumerate(ID):
		j = int(j)
		seenMutsThisGeneration.append(int(j))
		if (j not in MutationFrequenciesp1):      			#we haven't seen this mutation.
			MutationFrequenciesp1[j] = []
			MutationFrequenciesp2[j] = []
			MutationFrequenciesp3[j] = []
			MutationStartTimes[j] = [intervalCount,numAxis]
		MutationFrequenciesp1[j].append(float(MarineFreqs[i]))
		MutationFrequenciesp2[j].append(float(FreshwaterFreqs[i]))
		MutationFrequenciesp3[j].append(float(IntroducedFreqs[i]))
		
	for i in MutationFrequenciesp1:
		if ((i not in seenMutsThisGeneration) & (MutationStartTimes[i][1] == numAxis)):
			MutationStartTimes[i][1] = intervalCount
	
	intervalCount += 1

			
disDict1 = {k:MutationFrequenciesp1[k] for k in (2312765,7286,36974,419858,10933,22680,1553296,2327910,4638175)}
disDict2 = {k:MutationFrequenciesp2[k] for k in (2312765,7286,36974,419858,10933,22680,1553296,2327910,4638175)}
disDict3 = {k:MutationFrequenciesp3[k] for k in (2312765,7286,36974,419858,10933,22680,1553296,2327910,4638175)}

Traces = []

for i in disDict1:
	
	yAxis = [None] * numAxis
	yAxis[MutationStartTimes[i][0]:MutationStartTimes[i][1]] = MutationFrequenciesp1[i]	

#	if(len(yAxis) != numAxis):
#		print("Something went wrong")

	trace = go.Scatter(
	    x = xAxis,
	    y = yAxis,
	    name = str(i) + "p1",
	    line = dict(
		color = ("rgb(10,10,100)"),
		width = 1)
	)

	Traces.append(trace)
	
for i in disDict2:
	
	yAxis = [None] * numAxis
	yAxis[MutationStartTimes[i][0]:MutationStartTimes[i][1]] = MutationFrequenciesp2[i]	
	
	trace = go.Scatter(
	    x = xAxis,
	    y = yAxis,
	    name = str(i) + "p2",
	    line = dict(
		color = ("rgb(10,100,10)"),
		width = 1)
	)

	Traces.append(trace)
	

for i in disDict3:
	
	yAxis = [None] * numAxis
	yAxis[MutationStartTimes[i][0]:MutationStartTimes[i][1]] = MutationFrequenciesp3[i]	

	
	trace = go.Scatter(
	    x = xAxis,
	    y = yAxis,
	    name = str(i) + "p3",
	    line = dict(
		color = ("rgb(100,10,10)"),
		width = 1)
	)

	Traces.append(trace)

print(len(Traces))	

data = Traces




# Edit the layout
layout = dict(title = 'Effect Mutation Frequency Trajectory',
              xaxis = dict(title = 'Generations'),
              yaxis = dict(title = 'Frequency'),
              )

fig = dict(data=data, layout=layout)
#py.iplot(fig, filename='3_0_1 Average Fst')
plotly.offline.plot(fig, filename='8_0_8EffectMutationFrequencyTrajectory_FstAbove50.html')

