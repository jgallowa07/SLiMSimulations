initialize()
{
	
	numLakes = 10;
		
	defineConstant("NLAKES",numLakes);
	
	initializeSLiMOptions(dimensionality="x");
	
	initializeMutationRate(1e-7);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1",m1,1.0);
	initializeGenomicElement(g1, 0, 99999);
	initializeRecombinationRate(1e-8);
	
	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.1);
	i1.setInteractionFunction("n", 0.5, 0.02);
	
}

1 late(){
	sim.addSubpop("p1", 2000);
	IPF = 200;	
	p1.setSpatialBounds(c(0.0, 10.0));
	
	lakeIndividuals = p1.individuals;

	//Set 200 individuals in each 'freshwater lake'

	for (i in 0:(NLAKES-1)){
		for (j in 0:(IPF-1)){
			lakeIndividuals[IPF*i + j].setSpatialPosition((i+1)-(0.5));						
		}
	}
	
	catn(p1.individuals.spatialPosition);
}

1:late(){

	i1.evaluate();
		
	//sample how many individuals are in each lake
	
	catn("Generation: " +  sim.generation + "----------------------------------------------------------------");	
	for (l in 1:NLAKES){
		locationOfLake = (l - 0.5);
		lake_individuals = p1.individuals[abs(p1.individuals.spatialPosition - locationOfLake) < 1e-8];
	
		
		catn("Num at position "+locationOfLake+ ": " + size(lake_individuals) + "\n");	
		
		if(size(lake_individuals) == 0)
		{
			catn("DEADLAKEBUG?  :'("+ "\n");
			sim.simulationFinished();	
		}

		
		
	}
}
1: mateChoice() {
	
	//make sure second parent chosen is from the same 'freshwater lake'
	
	return i1.strength(individual); 
	//neighbors = i1.nearestNeighbors(individual, 3);
	//mates = rep(0.0, p1.individualCount);
	//mates[neighbors.index] = 1.0;
	//return mates;	
}

modifyChild(p1){

	//put child in the same lake as the parents
	 
	child.setSpatialPosition(parent1.spatialPosition);

	if (parent1.spatialPosition != child.spatialPosition | parent2.spatialPosition != child.spatialPosition)
	{
		catn("wtf" + "\n");
		sim.simulationFinished();
				

	} 

	return T;
}

1000 late ()
{

	sim.simulationFinished();
}
