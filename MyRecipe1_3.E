initialize() {
	initializeMutationRate(1e-5);
	
	scriptForQTLs = "if (runif(1) < 0.5) -1; else 1;";
	initializeMutationType("m1", 0.5, "s", scriptForQTLs);
	initializeGenomicElementType("g1", m1, 1.0);
	//m1.convertToSubstitution = F;
	//m1.mutationStackPolicy = "l";
	initializeGenomicElement(g1,0,999);
	initializeRecombinationRate(1e-3);
	end = 40001;
	interval = 40000;
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

		//cat(mean(p1.genomes.countOfMutationsOfType(m1))+ " " + mean(p2.genomes.countOfMutationsOfType(m1))+ "\n");	
		
		//cat(sim.mutations);	
		

		Chrom_Pos = 0:999;
		muts = sim.mutationsOfType(m1);
		polymorphic_sites = 0:0;
		//respective_frequenciesp1 = 0:0;
		//respective_frequenciesp2 = 0:0;
		numIndiv = size(sim.subpopulations.individuals);
		for (i in Chrom_Pos){
			number_of_muts_at_i = size(muts[muts.position == i]);
			if (number_of_muts_at_i > 0) {
				polymorphic_sites = c(polymorphic_sites,i);
				//respective_frequency = number_of_muts_at_i / numIndiv;
				//respective_frequencies = c(respective_frequencies,respective_frequency);
				
			}
		}
		//cat(size(polymorphic_sites)+" "+0+"\n");
		//cat(size(respective_frequencies)+"\n");
		//cat("\n");
		//cat("Size of Polymorphic_sites (after): " + size(polymorphic_sites) + "\n");
		//cat(polymorphic_sites);
		
		plus = muts[muts.selectionCoeff == 1.0];
		minus = muts[muts.selectionCoeff == -1.0];
		
		cat("\nOverall frequencies:\n\n");
		cat("Generation: "+ sim.generation + "\n");
		for (q in polymorphic_sites)
		{
			qPlus = plus[plus.position == q];
			qMinus = minus[minus.position == q];
			pf = sum(sim.mutationFrequencies(NULL, qPlus));
			mf = sum(sim.mutationFrequencies(NULL, qMinus));
			pf1 = sum(sim.mutationFrequencies(p1, qPlus));
			mf1 = sum(sim.mutationFrequencies(p1, qMinus));
			pf2 = sum(sim.mutationFrequencies(p2, qPlus));
			mf2 = sum(sim.mutationFrequencies(p2, qMinus));
			
			cat("   QTL " + q + ": f(+) == " + pf + ", f(-) == " + mf + "\n");
			cat("         in p1: f(+) == " + pf1 + ", f(-) == " + mf1 + "\n");
			cat("         in p2: f(+) == " + pf2 + ", f(-) == " + mf2 + "\n\n");
		} 
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

40001 late() {sim.simulationFinished();}

