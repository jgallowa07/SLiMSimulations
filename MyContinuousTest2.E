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

	for (i in 0:(NLAKES-1)){
		for (j in 0:(IPF-1)){
			lakeIndividuals[IPF*i + j].setSpatialPosition((i+1)-(0.5));						
		}
	}
	
	catn(p1.individuals.spatialPosition);
}

1:late(){

	i1.evaluate();

	catn("Generation: " +  sim.generation + "----------------------------------------------------------------");	
	for (l in 1:NLAKES){
		locationOfLake = (l - 0.5);
		lake_individuals = p1.individuals[abs(p1.individuals.spatialPosition - locationOfLake) < 1e-8];
	
		
		catn("Num at position "+locationOfLake+ ": " + size(lake_individuals) + "\n");	
		
		if(size(lake_individuals) == 0)
		{
			catn("DEADLAKE:'("+ "\n");
			sim.simulationFinished();	
		}

		
		
	}
}
1: mateChoice() {
	return i1.strength(individual); 
}

modifyChild(p1){
	child.setSpatialPosition(parent1.spatialPosition);
	return T;
}

1000 late ()
{

	sim.simulationFinished();
}
