initialize() {


	endOfSimulation = 1000;
	samplingInterval = 100;

	numberOfLakes = 10;	

	individualsInOcean = 2000;
	
	individualsInFresh = 2000;

	individualsPerFresh = individualsInFresh/numberOfLakes;	

	lengthOfPops = asFloat(numberOfLakes);
	numberOfSamplingPoints = endOfSimulation/samplingInterval + 1;
		
	defineConstant("END",endOfSimulation);
	defineConstant("INT",samplingInterval);
	defineConstant("NLAKES",numberOfLakes);
	defineConstant("IIO", individualsInOcean);
	defineConstant("IIF", individualsInFresh);	
	defineConstant("LOP", lengthOfPops);
	defineConstant("IPF",individualsPerFresh);

	initializeSLiMOptions(dimensionality="x");
	initializeMutationRate(1e-5);	
	initializeMutationType("m1", 0.5, "g", 0.5, 1.0);
	initializeMutationType("m2", 0.5, "g", -0.5, 1.0);
	initializeGenomicElementType("g1", c(m1,m2), c(1.0,1.0));
	m1.mutationStackPolicy = "l";
	m2.mutationStackPolicy = "l";
	initializeGenomicElement(g1,0,999);
	initializeRecombinationRate(1e-3);
	
	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.02);
	i1.setInteractionFunction("n", 0.5, 0.02);
	
}


1 late(){
	sim.addSubpop("p1", IIO);
	sim.addSubpop("p2", IIF);
	
	p1.setSpatialBounds(c(0.0, LOP));
	p2.setSpatialBounds(c(0.0, LOP));
	
	p1.setMigrationRates(p2,0.0005);
	p2.setMigrationRates(p1,0.0005);
	
	for (ind in p1.individuals){
		ind.setSpatialPosition(p1.pointUniform());
	}
	
	lakeIndividuals = p2.individuals;

	for (i in 0:(NLAKES-1)){
		for (j in 0:(IPF-1)){
			lakeIndividuals[IPF*i + j].setSpatialPosition((i+1)-(0.5));						
		}
	}

}
	

1 late() {

	//choose 100 locations along the chromosome, give each genome a mutation of type m1 or m2 at that location
 
	g = sim.subpopulations.genomes;
	n = size(g);
	mutations  = asInteger(round(runif(10,0,999)));
	mutationSites = sort(mutations);

	//for (m in mutationSites)
	//{
	//	isPlus = asLogical(rbinom(n, 1, 0.5));
	//	g[isPlus].addNewMutation(m1, 1.0, m);
	//	g[!isPlus].addNewMutation(m2, -1.0, m);
	//}
	
}


1: late() {
	i1.evaluate();
		
	Oceanic_mean_phenotype = NULL;
	Freshwater_mean_phenotype = NULL;
	substitutionOffset  = 0;

	for (fix in sim.substitutions){
		substitutionOffset = substitutionOffset + fix.selectionCoeff;

	}	

	for (individual in sim.subpopulations.individuals)
	{
		individual.tagF = substitutionOffset + sum(individual.genomes.mutations.selectionCoeff);		
	}
	

	//The entire if block below is for sampling data at a given interval, INT	
	
	if(((sim.generation-1) % INT == 0) | (sim.generation == END))
	{
		cat("-------------------------------\n");
		cat("Output for generation " + sim.generation + ":\n\n");
	}
	if(sim.generation == END)
	{;
		
		muts = sim.mutations;
		num = size(muts);
		mutfreq = sim.mutationFrequencies(NULL,muts);
		muteff = muts.selectionCoeff;
		mutpos = muts.position;
		writeFile("./myOut194_Muts.txt",num+" ",append = F);
		for (i in 0:(num-1))
		{
			writeFile("./myOut194_Muts.txt",mutfreq[i]+" "+mutpos[i]+" "+muteff[i],append = T);
		}
					


		sim.simulationFinished();
	}

	//Set tag equal to a value to be returned by fitness callback 
	
	for (individual in p2.individuals){
		individual.tagF = dnorm(-10.0 - individual.tagF,0.0,10.0);
	}
	for (individual in p1.individuals){
		individual.tagF = dnorm(10 - individual.tagF,0.0,10.0);
	}

	//Normalize the fitness in the lakes so an equivelent mean number of parents get chosen from each.
 
	for (l in 1:NLAKES){
		locationOfLake = (l - 0.5);
		lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
		sumOfPhenos = sum(lake_individuals.tagF);
		lake_individuals.tagF =(lake_individuals.tagF*IPF)/(sumOfPhenos);
	}
	
}


fitness(m1) {return 1.0;}
fitness(m2) {return 1.0;}
fitness(NULL, p1) {
	//phenotype = asFloat(individual.tagF);
	//return 1.0 + dnorm(10.0 - phenotype, 0.0, 10.0);    // optimum of +10
	return 1.0 + individual.tagF;
}
fitness(NULL, p2) {
	//phenotype = asFloat(individual.tagF);
	//return 1.0 + dnorm(- 10.0 - phenotype, 0.0, 10.0);    // optimum of +10
	return 1.0 + individual.tagF;
}

1: mateChoice() {
	return i1.strength(individual); 
}

modifyChild(p1){

	if(sourceSubpop != subpop){
		num = parent1.spatialPosition + rnorm(1,0,1.5);
		num = p2.pointReflected(num);
		//nearestMult = num + H - (num + H) % M;
		nearestMult = round(num);
		if(num < nearestMult){
			child.setSpatialPosition(nearestMult - 0.5);
		}else{
			child.setSpatialPosition(nearestMult + 0.5);		
		}
		
	}else{
		pos = parent1.spatialPosition + rnorm(1,0,1.5);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}
	return T;
}

modifyChild(p2){
	if(sourceSubpop != subpop){
		pos = parent1.spatialPosition + rnorm(1,0,0.5);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}else{
		child.setSpatialPosition(parent1.spatialPosition);
	}
	return T;
}


40001 late() {sim.simulationFinished();}
