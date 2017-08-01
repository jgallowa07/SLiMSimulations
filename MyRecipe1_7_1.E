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
	interval = 20000;
	numPoint = end/interval;
	defineConstant("END",end);
	defineConstant("Int",interval);
	defineConstant("numP",numPoint);	
	
}
1 early() {
	subPopSize = 10;
	sim.addSubpop("p1", 1800); 
	freshWater = 2:subPopSize;
	for (i in freshWater)
	{
		//print("iteration"+i);
		//f = sample(1:300,1);
		//print(f);
		sim.addSubpop(i,200);
	}
	subPops = sim.subpopulations;

	x = 1:(subPopSize-1);
	
	for (j in x)
	{
		subPops[j].setMigrationRates(1,0.0005);
	}
	for (k in freshWater)
	{
		subPops[0].setMigrationRates(k,0.005);
	}

	//sim.tag = 2;
	//sim.addSubpop("p1", 1000);
	//sim.addSubpop("p2", 1000);
	
	// set up migration; comment these out for zero gene flow
	//p1.setMigrationRates(p2, 0.01);
	//p2.setMigrationRates(p1, 0.01);
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
	
	high_freq_Plus = NULL;
	high_freq_Minus = NULL;
	low_freq_Plus = NULL;
	low_freq_Minus = NULL;
	
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
	
	if(((sim.generation-1) % Int == 0) | (sim.generation == END))
	{
		//cat("-------------------------------\n");
		//cat("Output for generation " + sim.generation + ":\n\n");
		p1_mean_fitness = mean(p1.cachedFitness(NULL));
		p2_mean_fitness = mean(p2.cachedFitness(NULL));
		p1_mean_phenotype =  mean(p1.individuals.tag);
		p2_mean_phenotype = mean(p2.individuals.tag);
	
		//cat(mean(p1.genomes.countOfMutationsOfType(m1))+ " " + mean(p2.genomes.countOfMutationsOfType(m1))+ "\n");	
		
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
		//cat("hard: "+size(polymorphic_sites)+"\n");
		//cat("easy: "+polymorphic_sites_easy+"\n");
		//cat(size(polymorphic_sites) + " " + 0 + "\n");	
		//cat(size(polymorphic_sites)+" "+p1_mean_fitness+" "+p2_mean_fitness+" "+p1_mean_phenotype+" "+p2_mean_phenotype+"\n");

		for (q in polymorphic_sites)
		{
			qPlus = m1_muts[m1_muts.position == q];
			qMinus = m2_muts[m2_muts.position == q];
			pf = sum(sim.mutationFrequencies(NULL, qPlus));
			mf = sum(sim.mutationFrequencies(NULL, qMinus));
		
			if (pf > 0.2 & pf < 0.8){ 
				high_freq_Plus = c(high_freq_Plus,q);
				qPlus.tag = 0;
			}else{
				low_freq_Plus = c(low_freq_Plus,q);
				qPlus.tag = 1;
			}
			if (mf > 0.2 & mf < 0.8){ 
				high_freq_Minus = c(high_freq_Minus,q);
				qMinus.tag = 2;						
			}else{
				low_freq_Minus = c(low_freq_Minus,q);
				qMinus.tag = 3;
			}
		}
		numHFP = size(high_freq_Plus);
		numLFP = size(low_freq_Plus);
		numHFM = size(high_freq_Minus);
		numLFM = size(low_freq_Minus);
	}
	if(sim.generation == END)
	{

		individuals = sim.subpopulations.individuals;
		
		respectiveNumHighFreq_m1 = NULL;			
		respectiveNumLowFreq_m1 = NULL;			
		respectiveNumHighFreq_m2 = NULL;			
		respectiveNumLowFreq_m2 = NULL;			
		
		
	
		for (i in individuals){
			
			m1mutations = i.uniqueMutationsOfType(m1);
			m2mutations = i.uniqueMutationsOfType(m2);

			indNumOfHighFreq_m1 = 0;
			indNumOfLowFreq_m1 = 0;
			indNumOfHighFreq_m2 = 0;
			indNumOfLowFreq_m2 = 0;
			
			for (m_1 in m1mutations){
				//Funky logic alert y'alls
				//if true:
				//this means that the m_1 mutation position is contained in the vector of high_freq_Plus (plus == m1)
				//i.e m_1 is a high frequency mutation. 
				if(size(high_freq_Plus[high_freq_Plus == m_1.position]) > 0){	
					indNumOfHighFreq_m1 = indNumOfHighFreq_m1 + 1;
				}else{
					indNumOfLowFreq_m1 = indNumOfLowFreq_m1 + 1;
				}
			}	
			for (m_2 in m2mutations){	
				if(size(high_freq_Minus[high_freq_Minus == m_2.position]) > 0){	
					indNumOfHighFreq_m2 = indNumOfHighFreq_m2 + 1;
				}else{
					indNumOfLowFreq_m2 = indNumOfLowFreq_m2 + 1;
				}
			}
			respectiveNumHighFreq_m1 = c(respectiveNumHighFreq_m1,indNumOfHighFreq_m1);
			respectiveNumLowFreq_m1 = c(respectiveNumLowFreq_m1,indNumOfLowFreq_m1);
			respectiveNumHighFreq_m2 = c(respectiveNumHighFreq_m2,indNumOfHighFreq_m2);
			respectiveNumLowFreq_m2 = c(respectiveNumLowFreq_m2,indNumOfLowFreq_m2);
			
			//catn("OLDWAY--------------------------------");
			cat(indNumOfHighFreq_m1+" "+indNumOfLowFreq_m1+" "+indNumOfHighFreq_m2+" "+indNumOfLowFreq_m2+"\n");
		}
		//catn("num high freq orig: " + size(respectiveNumHighFreq_m1));
		//nhf = sim.mutations[sim.mutations.tag == 0]; 
		//catn("num high freq new: "+ size(nhf));	
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
