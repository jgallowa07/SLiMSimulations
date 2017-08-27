initialize() {
	initializeMutationRate(1e-5);
	
	initializeMutationType("m1", 0.5, "f", 1.0);
	initializeMutationType("m2", 0.5, "f", -1.0);
	initializeGenomicElementType("g1", c(m1,m2), c(1.0,1.0));
	m1.mutationStackPolicy = "l";
	m2.mutationStackPolicy = "l";
	initializeGenomicElement(g1,0,999);
	initializeRecombinationRate(1e-3);

	end = 20000;
	interval = 50;
	numPoint = end/interval;
	subPopSize = 11;
	//freshWater = subPopSize - 1;
	defineConstant("SubPopSize", subPopSize);
	defineConstant("END",end);
	defineConstant("Int",interval);
	defineConstant("numP",numPoint);	
	
	
}
1 early() {
	//subPopSize = SubPopSize;
	sim.addSubpop("p1", 2000); 
	freshWater = 2:SubPopSize;
	for (i in freshWater)
	{
		//f = sample(1:300,1);
		//print(f);
		sim.addSubpop(i,200);
	}
	subPops = sim.subpopulations;

	x = 1:(SubPopSize-1);
	
	for (j in x)
	{
		subPops[j].setMigrationRates(1,0.0005);
	}
	for (k in freshWater)
	{
		subPops[0].setMigrationRates(k,0.00005);
	}
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
				//a = rexp(1,0.5);
				//b = rexp(5);
				//x = a / (a+b);
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
		
	Oceanic_mean_phenotype = NULL;
	Freshwater_mean_phenotype = NULL;
	
	
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
		Oceanic_mean_phenotype =  mean(p1.individuals.tag);
		mean_fresh = NULL;
		lastfreshindex = SubPopSize - 1;
		for (i in 1:lastfreshindex){
		
			mean_fresh = c(mean_fresh,mean(sim.subpopulations[i].individuals.tag));
		}
		Freshwater_mean_phenotype = mean(mean_fresh);

		
		writeFile("./myOut175_Ph.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype,append = T);
	
		//cat(mean(p1.genomes.countOfMutationsOfType(m1))+ " " + mean(p2.genomes.countOfMutationsOfType(m1))+ "\n");	
		

			
		
	}
	if(sim.generation == END)
	{

		Chrom_Pos = 0:999;
		m1_muts = sim.mutationsOfType(m1);
		m2_muts = sim.mutationsOfType(m2);
		polymorphic_sites = NULL;
	
		//polymorphic_sites_easy = unique(sim.mutations.position);		

		for (i in Chrom_Pos){
			number_of_m1muts_at_i = size(m1_muts[m1_muts.position == i]);
			number_of_m2muts_at_i = size(m2_muts[m2_muts.position == i]);
			if (number_of_m1muts_at_i > 0 | number_of_m2muts_at_i > 0) {
				polymorphic_sites = c(polymorphic_sites,i);
			}
		}
		//cat(size(polymorphic_sites) + " " + 0 + "\n");	
		//cat(size(polymorphic_sites)+" "+p1_mean_fitness+" "+p2_mean_fitness+" "+p1_mean_phenotype+" "+p2_mean_phenotype+"\n");

		for (q in polymorphic_sites)
		{
			qPlus = m1_muts[m1_muts.position == q];
			qMinus = m2_muts[m2_muts.position == q];
			pf = sum(sim.mutationFrequencies(NULL, qPlus));
			mf = sum(sim.mutationFrequencies(NULL, qMinus));

			if (pf > 0.2 & pf < 0.8){ 
				qPlus.tag = 0;
			}else{
				qPlus.tag = 1;
			}
			if (mf > 0.2 & mf < 0.8){ 
				qMinus.tag = 2;						
			}else{
				qMinus.tag = 3;
			}
		}
		Oceanic = p1.individuals;
		lastfreshindex = SubPopSize - 1;
		Freshwater = 1:lastfreshindex;			
	
		for (n in Oceanic){
			
			m1mutations = n.uniqueMutationsOfType(m1);
			m2mutations = n.uniqueMutationsOfType(m2);
			//Oinoh/lm1/2 = "oceanic individual's # of high/low frequency m1/m2"
			Oinohfm1 = size(m1mutations[m1mutations.tag == 0]);
			Oinolfm1 = size(m1mutations[m1mutations.tag == 1]);
	
			Oinohfm2 = size(m2mutations[m2mutations.tag == 2]); 	
			Oinolfm2 = size(m2mutations[m2mutations.tag == 3]);

			writeFile("./myOut175_Oc.txt",Oinohfm1+" "+Oinolfm1+" "+Oinohfm2+" "+Oinolfm2,append = T);
		
		} 

		for (subPop in Freshwater){

			Finohfm1 = NULL;
			Finolfm1 = NULL;
			Finohfm2 = NULL;
			Finolfm2 = NULL;
		
			for (individual in sim.subpopulations[subPop].individuals){
				m1muts = individual.uniqueMutationsOfType(m1);
				m2muts = individual.uniqueMutationsOfType(m2);				

				Finohfm1 = c(Finohfm1,size(m1muts[m1muts.tag == 0]));
				Finolfm1 = c(Finolfm1,size(m1muts[m1muts.tag == 1]));
				Finohfm2 = c(Finohfm2,size(m2muts[m2muts.tag == 2]));
				Finolfm2 = c(Finolfm2,size(m2muts[m2muts.tag == 3]));
			}
			
			for (z in 0:199){
				writeFile("./myOut175_Fr.txt",Finohfm1[z]+" "+Finolfm1[z]+" "+Finohfm2[z]+" "+Finolfm2[z],append = T);
			}			
		}			
	
		//catn("Oceanic mean phenotype: " + Oceanic_mean_phenotype);
		//catn("Freshwater Mean Phenotype: " + Freshwater_mean_phenotype);
		
		sim.simulationFinished();
	}
}


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
fitness(NULL, p3) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}
fitness(NULL, p4) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}
fitness(NULL, p5) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}
fitness(NULL, p6) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}
fitness(NULL, p7) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}
fitness(NULL, p8) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}
fitness(NULL, p9) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}
fitness(NULL, p10) {
	phenotype = asFloat(individual.tag);
	return 1.0 + dnorm(-10.0 - phenotype, 0.0, 5.0);   // optimum of -10
}
40001 late() {sim.simulationFinished();}
