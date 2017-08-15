initialize() {


	end = 20000;
	interval = 250;

	numLakes = 10;	

	IndividualsInOcean = 2000;
	
	IndividualsInFresh = 200 * numLakes;	

	lengthOfPops = asFloat(numLakes);
	numPoint = end/interval;
		

	defineConstant("END",end);
	defineConstant("INT",interval);
	defineConstant("NLAKES",numLakes);
	defineConstant("IIO", IndividualsInOcean);
	defineConstant("IIF", IndividualsInFresh);	
	defineConstant("LOP", lengthOfPops);
	defineConstant("IPF",IndividualsInFresh/numLakes);
	initializeSLiMOptions(dimensionality="x");

	initializeMutationRate(1e-5);
	
	initializeMutationType("m1", 0.5, "f", 1.0);
	initializeMutationType("m2", 0.5, "f", -1.0);
	initializeGenomicElementType("g1", c(m1,m2), c(1.0,1.0));
	m1.mutationStackPolicy = "l";
	m2.mutationStackPolicy = "l";
	initializeGenomicElement(g1,0,999);
	initializeRecombinationRate(1e-3);
	
	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.02);
	i1.setInteractionFunction("n", 0.5, 0.02);
	
	//initializeInteractionType(2, "x", reciprocal=T, maxDistance=0.3);
	//i2.setInteractionFunction("n",1.0,0.1);
	
	writeFile("./myOut182_Ph.txt",numPoint+" ",append = F);
	writeFile("./myOut182_Oc.txt",numPoint+" ",append = F);
	writeFile("./myOut182_Fr.txt",numPoint+" ",append = F);
	writeFile("./myOut182_Wh.txt",numPoint+" ",append = F);
}

//2: early() {	
//	for (l in 1:NLAKES){
//		locationOfLake = (l - 0.5);
//		lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
//		sumOfPhenos = sum(lake_individuals.tagF);
//		lake_individuals.tagF =(lake_individuals.tagF*IPF)/(sumOfPhenos + 1);
//	}
//}

1 late(){
	sim.addSubpop("p1", IIO);
	sim.addSubpop("p2", IIF);
	
	//IPF = IIF/NLAKES;

	p1.setSpatialBounds(c(0.0, LOP));
	p2.setSpatialBounds(c(0.0, LOP));
	
	p1.setMigrationRates(p2,0.0005);
	p2.setMigrationRates(p1,0.0005);
	
	for (ind in p1.individuals){
		ind.setSpatialPosition(p1.pointUniform());
	}
	
	lakeIndividuals = p2.individuals;

	//for (j in 0:(IPF-1)){
	for (i in 0:(NLAKES-1)){
		for (j in 0:(IPF-1)){
		//for (i in 0:(NLAKES-1)){
			lakeIndividuals[IPF*i + j].setSpatialPosition((i+1)-(0.5));						
		}
	}
	
	//catn(p2.individuals.spatialPosition);
}
//1: late(){
	//i1.evaluate();
	//catn("Generation: " +  sim.generation);	
	//Loc = 0.5;
	//for(i in 1:NLAKES){
	//	x = sum(abs(p2.individuals.spatialPosition - Loc) < 1e-8);
	//	catn("Num at position "+Loc+ ": " + x + "\n");	
	//	Loc = Loc + 1;
	//}
	

//}

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
	i1.evaluate();
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

	for (individual in p1.individuals)
	{
		m1muts = individual.countOfMutationsOfType(m1);
		m2muts = individual.countOfMutationsOfType(m2);	
		individual.tagF = substitutionOffset + (m1muts - m2muts);
	}
	for (individual in p2.individuals){
		m1muts = individual.countOfMutationsOfType(m1);
		m2muts = individual.countOfMutationsOfType(m2);	
		individual.tagF = dnorm(-10.0 - substitutionOffset + (m1muts - m2muts),0.0,10.0);
	}

	//The entire if block below is for sampling data at a given interval, INT	
	
	if(((sim.generation-1) % INT == 0) | (sim.generation == END))
	{
		cat("-------------------------------\n");
		cat("Output for generation " + sim.generation + ":\n\n");
		p1_mean_fitness = mean(p1.cachedFitness(NULL));
		p2_mean_fitness = mean(p2.cachedFitness(NULL));
		Oceanic_mean_phenotype =  mean(p1.individuals.tagF);
		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
	
		writeFile("./myOut182_Ph.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype,append = T);
		
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
		Freshwater = p2.individuals;

		oceanicMeanNumberOfHighFreqM1 = NULL;
		oceanicMeanNumberOfLowFreqM1 = NULL;
		oceanicMeanNumberOfHighFreqM2 = NULL;
		oceanicMeanNumberOfLowFreqM2 = NULL;

		freshwaterMeanNumberOfHighFreqM1 = NULL;
		freshwaterMeanNumberOfLowFreqM1 = NULL;
		freshwaterMeanNumberOfHighFreqM2 = NULL;
		freshwaterMeanNumberOfLowFreqM2 = NULL;
	

		for (o_individual in Oceanic){
			
			muts = o_individual.genomes.mutations;

			//Oinoh/lm1/2 = "oceanic individual's # of high/low frequency m1/m2"
			oceanicMeanNumberOfHighFreqM1 = c(oceanicMeanNumberOfHighFreqM1,size(muts[muts.tag == 0]));
			oceanicMeanNumberOfLowFreqM1 = c(oceanicMeanNumberOfLowFreqM1,size(muts[muts.tag == 1]));
			oceanicMeanNumberOfHighFreqM2 = c(oceanicMeanNumberOfHighFreqM2,size(muts[muts.tag == 2])); 	
			oceanicMeanNumberOfLowFreqM2 = c(oceanicMeanNumberOfLowFreqM2,size(muts[muts.tag == 3]));
		}
 
		for (f_individual in Freshwater){

			mutss = f_individual.genomes.mutations;

			freshwaterMeanNumberOfHighFreqM1 = c(freshwaterMeanNumberOfHighFreqM1,size(mutss[mutss.tag == 0]));
			freshwaterMeanNumberOfLowFreqM1 = c(freshwaterMeanNumberOfLowFreqM1,size(mutss[mutss.tag == 1]));
			freshwaterMeanNumberOfHighFreqM2 = c(freshwaterMeanNumberOfHighFreqM2,size(mutss[mutss.tag == 2]));
			freshwaterMeanNumberOfLowFreqM2 = c(freshwaterMeanNumberOfLowFreqM2,size(mutss[mutss.tag == 3]));

		}
	
		omnohfm1 = mean(oceanicMeanNumberOfHighFreqM1);
                omnolfm1 = mean(oceanicMeanNumberOfLowFreqM1);
                omnohfm2 = mean(oceanicMeanNumberOfHighFreqM2);
                omnolfm2 = mean(oceanicMeanNumberOfLowFreqM2);

                fmnohfm1 = mean(freshwaterMeanNumberOfHighFreqM1);
                fmnolfm1 = mean(freshwaterMeanNumberOfLowFreqM1);
                fmnohfm2 = mean(freshwaterMeanNumberOfHighFreqM2);
                fmnolfm2 = mean(freshwaterMeanNumberOfLowFreqM2);

		wmnohfm1 = mean(c(omnohfm1,fmnohfm1));
		wmnolfm1 = mean(c(omnolfm1,fmnolfm1));
		wmnohfm2 = mean(c(omnohfm2,fmnohfm2));
		wmnolfm2 = mean(c(omnolfm2,fmnolfm2));

		writeFile("./myOut182_Wh.txt",wmnohfm1+" "+wmnolfm1+" "+wmnohfm2+" "+wmnolfm2,append = T);
		writeFile("./myOut182_Oc.txt",omnohfm1+" "+omnolfm1+" "+omnohfm2+" "+omnolfm2,append = T);
		writeFile("./myOut182_Fr.txt",fmnohfm1+" "+fmnolfm1+" "+fmnohfm2+" "+fmnolfm2,append = T);
			
	}
	if(sim.generation == END)
	{
		sim.simulationFinished();
//		catn(p2.individuals.spatialPosition);
	}
//	catn("Generation: " +  sim.generation + "----------------------------------------------------------------");	
	for (l in 1:NLAKES){
		locationOfLake = (l - 0.5);
		lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
		sumOfPhenos = sum(lake_individuals.tagF);
		lake_individuals.tagF =(lake_individuals.tagF*IPF)/(sumOfPhenos);
	
		
//		catn("Num at position "+locationOfLake+ ": " + size(lake_individuals) + "\n");	
//		x =NULL;
//		if (size(lake_individuals) == 0){
//			x = 0;
//		}else{
//			x = mean(lake_individuals.tagF);
//		}  
//		catn("Average Phenotype (fitness) at location " + locationOfLake + ":  " + x + "\n");
		
	}
	
}


fitness(m1) {return 1.0;}
fitness(m2) {return 1.0;}
fitness(NULL, p1) {
	phenotype = asFloat(individual.tagF);
	return 1.0 + dnorm(10.0 - phenotype, 0.0, 10.0);    // optimum of +10
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
		pos = parent1.spatialPosition + rnorm(1,0,0.5);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}
	return T;
}

modifyChild(p2){
	if(sourceSubpop != subpop){
		pos = parent1.spatialPosition + rnorm(1,0,1.5);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}else{
		child.setSpatialPosition(parent1.spatialPosition);
	}
	return T;
}


40001 late() {sim.simulationFinished();}
