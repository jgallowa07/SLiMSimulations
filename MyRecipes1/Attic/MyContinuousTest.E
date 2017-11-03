initialize()
{
	
	numLakes = 10;
	
	multiples = 1/numLakes;
	half = multiples / 2;
	
	defineConstant("NUMLakes",numLakes);
	defineConstant("M", multiples);
	defineConstant("H", half);
	
	initializeSLiMOptions(dimensionality="x");
	
	initializeMutationRate(1e-7);
	initializeMutationType("m1", 0.5, "f", 0.0);
	initializeGenomicElementType("g1",m1,1.0);
	initializeGenomicElement(g1, 0, 99999);
	initializeRecombinationRate(1e-8);
	
	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.1);
	i1.setInteractionFunction("n", 0.5, 0.02);
	
}

1 early(){
	sim.addSubpop("p1", 2000);
	sim.addSubpop("p2", 2000);
	
	p1.setSpatialBounds(c(0.0, 1.0));
	p2.setSpatialBounds(c(0.0, 1.0));
	
	p1.setMigrationRates(p2,0.0005);
	p2.setMigrationRates(p1,0.0005);
	
	for (ind in p1.individuals){
		ind.setSpatialPosition(p1.pointUniform());
	}
	
	lakeIndividuals = p2.individuals;

	for (i in 0:(NUMLakes-1)){
		for (j in 0:199){
			lakeIndividuals[200*i + j].setSpatialPosition(asFloat((i/10)+0.05));						
		}
	}
	
	//x = size(p2.individuals.spatialPosition[p2.individuals.spatialPosition == 0.15]);
	//catn("Num at position 0.15: " + x);	
	catn(p2.individuals.spatialPosition);
}

1: early(){
//	catn("\n");
//	catn("Generation: "  + sim.generation);
//	catn(p2.individuals.spatialPosition);
	i1.evaluate();
}

1: mateChoice() {
// spatial mate choice
	return i1.strength(individual); 
}

modifyChild(p1){
	if(sourceSubpop != subpop){
		num = parent1.spatialPosition + rnorm(1,0,0.02);
		nearestMult = num + H - (num + H) % M;
		if(num < nearestMult){
			child.setSpatialPosition(nearestMult - H);
		}else{
			child.setSpatialPosition(nearestMult + H);		
		}
		
	}else{
		pos = parent1.spatialPosition + rnorm(1,0,0.02);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}
	return T;
}

modifyChild(p2){
	child.setSpatialPosition(parent1.spatialPosition);
	return T;
}


100 late()
{
	location = 0.05;
	for(i in 1:10) {
		catn(location);
		x1 = size(p2.individuals[p2.individuals.spatialPosition == location]);

		catn("Num at position" +location +"(after sim): " + x1);
		
		location = asFloat(location + 0.1);
	}
//	x2 = size(p2.individuals[p2.individuals.spatialPosition == 0.55]);	
//	catn(x2);
	cat("\n");
	catn(p2.individuals.spatialPosition);
	sim.simulationFinished();
}
