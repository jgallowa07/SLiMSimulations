initialize()
{
	
    setSeed(123);
	numLakes = 3;
		
	defineConstant("NLAKES",numLakes);
	
	initializeSLiMOptions(dimensionality="x");
	
	initializeMutationRate(0.0);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1",m1,1.0);
	initializeGenomicElement(g1, 0, 10);
	initializeRecombinationRate(0.0);
	
	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.1);
	i1.setInteractionFunction("n", 0.5, 0.5);
	
}

1 late(){
	IPF = 20;	
	sim.addSubpop("p1", NLAKES*IPF);
	p1.setSpatialBounds(c(0.0, NLAKES*1.0));
	
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
	
	catn("Generation: " +  sim.generation + " total indivs:" + size(p1.individuals) + "----------------------------------------");	
    catn(sort(p1.individuals.spatialPosition));
	for (l in 1:NLAKES){
		locationOfLake = (l - 0.5);
		lake_individuals = p1.individuals[abs(p1.individuals.spatialPosition - locationOfLake) < 1e-2];
	
		
		catn("Num at position "+locationOfLake+ ": " + size(lake_individuals) + "\n");	
		
		if(size(lake_individuals) == 0)
		{
			catn("DEADLAKEBUG?  :'("+ "\n");
            catn("Generation: " +  sim.generation + " total indivs:" + size(p1.individuals) + "------------------------------------");	
            catn(sort(p1.individuals.spatialPosition));
			sim.simulationFinished();	
		}

		
		
	}
}
1: mateChoice() {
	
	//make sure second parent chosen is from the same 'freshwater lake'
    out = i1.strength(individual);
    catn("Interactions of indiv at " + individual.spatialPosition);
    catn(out);
	
    // This casues the bug:
	return out;

    // This avoids the bug:
    // return 1.0 + out * 0.0;
}

modifyChild(p1){

	//put child in the same lake as the parents
	
	child.setSpatialPosition(parent1.spatialPosition);
	return T;
}

1000 late ()
{

	sim.simulationFinished();
}
