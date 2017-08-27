initialize() {
	initializeMutationRate(1e-5);
	
	initializeMutationType("m1", 0.5, "f", 1.0);
	initializeMutationType("m2", 0.5, "f", -1.0);
	initializeGenomicElementType("g1", c(m1,m2), c(1.0,1.0));
	//m1.convertToSubstitution = F;
	//m2.convertToSubstitution = F;
	m1.mutationStackPolicy = "l";
	m2.mutationStackPolicy = "l";
	initializeGenomicElement(g1,0,999);
	initializeRecombinationRate(1e-3);

	end = 20000;
	interval = 50;
	numPoint = end/interval;
	defineConstant("END",end);
	defineConstant("Int",interval);
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
					g.addNewMutation(m2,-1.0,i);
				}
			}
		}
	}	
}
1: late() {
	m1subs = 0;
	m2subs = 0;
	
	for (fix in sim.substitutions){
		if(fix.mutationType  == m1){
			m1subs = m1subs + 1;
		}else{
			m2subs = m2subs + 1;
		} 
	}	
	substitutionOffset  = m1subs - m2subs;

	for (individual in sim.subpopulations.individuals)
	{
		m1muts = individual.countOfMutationsOfType(m1);
		m2muts = individual.countOfMutationsOfType(m2);	
		individual.tag = substitutionOffset + (m1muts - m2muts);
	}
	
	if((sim.generation-1) % Int == 0)
	{
		//cat("-------------------------------\n");
		//cat("Output for generation " + sim.generation + ":\n\n");
		p1_mean_fitness = mean(p1.cachedFitness(NULL));
		p2_mean_fitness = mean(p2.cachedFitness(NULL));
		p1_mean_phenotype =  mean(p1.individuals.tag);
		p2_mean_phenotype = mean(p2.individuals.tag);
	
		//cat("p1 mean fitness = " + mean(p1.cachedFitness(NULL)) + "\n");
		//cat("p2 mean fitness = " + mean(p2.cachedFitness(NULL)) + "\n");
	
		//cat("p1 mean phenotype = " + mean(p1.individuals.tag) + "\n");
		//cat("p2 mean phenotype = " + mean(p2.individuals.tag) + "\n");
			
		//cat(mean(p1.individuals.tag)+" "+mean(p2.individuals.tag)+"\n");

		//cat(mean(p1.genomes.countOfMutationsOfType(m1))+ " " + mean(p2.genomes.countOfMutationsOfType(m1))+ "\n");	
		
		
		Chrom_Pos = 0:999;
		m1_muts = sim.mutationsOfType(m1);
		m2_muts = sim.mutationsOfType(m2);
		polymorphic_sites = NULL;
		for (i in Chrom_Pos){
			number_of_m1muts_at_i = size(m1_muts[m1_muts.position == i]);
			number_of_m2muts_at_i = size(m2_muts[m2_muts.position == i]);
			if (number_of_m1muts_at_i > 0 | number_of_m2muts_at_i > 0) {
				polymorphic_sites = c(polymorphic_sites,i);
			}
		}
		
		//cat(size(polymorphic_sites) + " " + 0 + "\n");	
		//cat(size(polymorphic_sites)+" "+p1_mean_fitness+" "+p2_mean_fitness+" "+p1_mean_phenotype+" "+p2_mean_phenotype+"\n");

		//cat("\nOverall frequencies:\n\n");
		//cat("Generation: "+ sim.generation + "\n");
		high_freq_Plus = NULL;
		high_freq_Minus = NULL;
		low_freq_Plus = NULL;
		low_freq_Minus = NULL;
		for (q in polymorphic_sites)
		{
			qPlus = m1_muts[m1_muts.position == q];
			qMinus = m2_muts[m2_muts.position == q];
			pf = sum(sim.mutationFrequencies(NULL, qPlus));
			mf = sum(sim.mutationFrequencies(NULL, qMinus));
		
			if (pf > 0.2 & pf < 0.8){ 
				high_freq_Plus = c(high_freq_Plus,q);
			}else{
				low_freq_Plus = c(low_freq_Plus,q);
			}
			if (mf > 0.2 & mf < 0.8){ 
				high_freq_Minus = c(high_freq_Minus,q);
			}else{
				low_freq_Minus = c(low_freq_Minus,q);
			}
			
			
			//pf1 = sum(sim.mutationFrequencies(p1, qPlus));
			//mf1 = sum(sim.mutationFrequencies(p1, qMinus));
			//pf2 = sum(sim.mutationFrequencies(p2, qPlus));
			///mf2 = sum(sim.mutationFrequencies(p2, qMinus));
			
			//cat("   QTL " + q + ": f(+) == " + pf + ", f(-) == " + mf + "\n");
			//cat("         in p1: f(+) == " + pf1 + ", f(-) == " + mf1 + "\n");
			//cat("         in p2: f(+) == " + pf2 + ", f(-) == " + mf2 + "\n\n");
			
		}
		numHFP = size(high_freq_Plus);
		numLFP = size(low_freq_Plus);
		numHFM = size(high_freq_Minus);
		numLFM = size(low_freq_Minus);

		cat(numHFP+" "+numLFP+" "+numHFM+" "+numLFM+"\n");
			
		//cat("num high freq plus: " + size(high_freq_Plus) + "\n"); 
		//cat("num low freq plus: " + size(low_freq_Plus) + "\n"); 
		//cat("num high freq minus: " + size(high_freq_Minus) + "\n"); 
		//cat("num low freq minus: " + size(low_freq_Minus) + "\n"); 
	}
	if(sim.generation == END)
	{
		sim.simulationFinished();
		//cat("-------------------------------\n");
		//cat("Output for generation " + sim.generation + ":\n\n");
	
		//cat("p1 mean fitness = " + mean(p1.cachedFitness(NULL)) + "\n");
		//cat("p2 mean fitness = " + mean(p2.cachedFitness(NULL)) + "\n");
	
		//cat("p1 mean phenotype = " + mean(p1.individuals.tag) + "\n");
		//cat("p2 mean phenotype = " + mean(p2.individuals.tag) + "\n");
	}
}

//200{p1.outputSample(10);}

fitness(m1) {return 1.0;}
fitness(m2) {return 1.0;}
fitness(NULL, p1) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(10.0 - phenotype, 0.0, 5.0);    // optimum of +10
}
fitness(NULL, p2) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}
40001 late() {sim.simulationFinished();}
