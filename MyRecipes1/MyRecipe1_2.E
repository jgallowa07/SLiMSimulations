initialize() {
	initializeMutationRate(1e-5);
	
	scriptForQTLs = "if (runif(1) < 0.5) -1; else 1;";
	initializeMutationType("m1", 0.5, "s", scriptForQTLs);
	initializeGenomicElementType("g1", m1, 1.0);
	//m1.convertToSubstitution = F;
	//m1.mutationStackPolicy = "l";
	initializeGenomicElement(g1,0,999);
	initializeRecombinationRate(1e-3);
	end = 40000;
	interval = 10;
	numPoint = end/interval;
	defineConstant("END",end);
	defineConstant("INT",interval);
	defineConstant("numP",numPoint);	
	//cat(numP+"\n");	
	
}
1 early() {
	sim.tag = 2;
	sim.addSubpop("p1", 1000);
	sim.addSubpop("p2", 1000);
	
	// set up migration; comment these out for zero gene flow
	p1.setMigrationRates(p2, 0.01);
	p2.setMigrationRates(p1, 0.01);
	cat("DATA:"+"\n");
	cat(numP+"\n");	
	
}

1 late() {
	// optional: give m1 mutations to everyone, as standing variation
	// if this is commented out, QTLs effectively start as 0
	G = sim.subpopulations.genomes;
	n = size(G);
	X = 0:999;
	
	Mutations = asLogical(rbinom(1000,1,0.1));
	for (i in X){
		if(Mutations[i]){
			for (g in G){	
				if(runif(1) < 0.5){
					g.addNewMutation(m1,1.0,i);
				}else{
					g.addNewMutation(m1,-1.0,i);
				}
			}
		}
	}	
			
				
}
1: late() {
	// construct phenotypes for the additive effects of QTLs
	for (individual in sim.subpopulations.individuals)
	{
		m1muts = individual.genomes.mutationsOfType(m1);
		individual.tag = asInteger(sum(m1muts.selectionCoeff));
	}
	
	if((sim.generation-1) % INT == 0)
	{
			//cat("-------------------------------\n");
			//cat("Output for generation " + sim.generation + ":\n\n");
	
			//cat("p1 mean fitness = " + mean(p1.cachedFitness(NULL)) + "\n");
			//cat("p2 mean fitness = " + mean(p2.cachedFitness(NULL)) + "\n");
	
			//cat("p1 mean phenotype = " + mean(p1.individuals.tag) + "\n");
			//cat("p2 mean phenotype = " + mean(p2.individuals.tag) + "\n");
			
		//cat(mean(p1.individuals.tag)+" "+mean(p2.individuals.tag)+"\n");

		cat(mean(p1.genomes.countOfMutationsOfType(m1))+ " " + mean(p2.genomes.countOfMutationsOfType(m1))+ "\n");	
			
	}
	if(sim.generation == END)
	{
		sim.simulationFinished();
	}
	
}

//200{p1.outputSample(10);}

fitness(m1) {
	// the QTLs themselves are neutral; their effect is handled below
	return 1.0;
}
fitness(NULL, p1) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(10.0 - phenotype, 0.0, 5.0);    // optimum of +10
}
fitness(NULL, p2) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}

40000 late() {sim.simulationFinished();}

