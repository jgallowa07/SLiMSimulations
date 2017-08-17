initialize() {


	endOfSimulation = 200;
	samplingInterval = 50;

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
	
	writeFile("./myOut193_Ph.txt",numberOfSamplingPoints+" "+samplingInterval,append = F);
	writeFile("./myOut193_Oc.txt",numberOfSamplingPoints+" "+samplingInterval,append = F);
	writeFile("./myOut193_Fr.txt",numberOfSamplingPoints+" "+samplingInterval,append = F);
	writeFile("./myOut193_Wh.txt",numberOfSamplingPoints+" "+samplingInterval,append = F);
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

	//catn(p2.individuals.spatialPosition);
}
	

1 late() {

	//choose 100 locations along the chromosome, give each genome a mutation of type m1 or m2 at that location
 
	g = sim.subpopulations.genomes;
	n = size(g);
	mutations  = asInteger(round(runif(100,0,999)));
	mutationSites = sort(mutations);

	for (m in mutationSites)
	{
		isPlus = asLogical(rbinom(n, 1, 0.5));
		g[isPlus].addNewMutation(m1, 1.0, m);
		g[!isPlus].addNewMutation(m2, -1.0, m);
	}
	
	//catn(sim.mutations.position);		
	//cat('\n');
	//catn(unique(sim.mutations.position));
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
		//p1_mean_fitness = mean(p1.cachedFitness(NULL));
		//p2_mean_fitness = mean(p2.cachedFitness(NULL));
		Oceanic_mean_phenotype =  mean(p1.individuals.tagF);
		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
		catn("SubstitutionOffset: "+substitutionOffset);
		catn("Oceanic Mean Phenotype: "+ Oceanic_mean_phenotype);
		catn("Freshwater Mean Phenotype: "+ Freshwater_mean_phenotype+"\n");
		writeFile("./myOut193_Ph.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype,append = T);
		
		Chrom_Pos = 0:999;
		m1_muts = sim.mutationsOfType(m1);
		m2_muts = sim.mutationsOfType(m2);
		//polymorphic_sites = NULL;
			
		//catn(size(sim.mutations));		
						
		polymorphic_sites = unique(sim.mutations.position);		

		//for (i in Chrom_Pos){
		//	number_of_m1muts_at_i = size(m1_muts[m1_muts.position == i]);
		//	number_of_m2muts_at_i = size(m2_muts[m2_muts.position == i]);
		//	if (number_of_m1muts_at_i > 0 | number_of_m2muts_at_i > 0) {
		//		polymorphic_sites = c(polymorphic_sites,i);
		//	}
		//}

		//catn("Hard way: "+size(polymorphic_sites));		
		//catn("Easy way: "+size(polymorphic_sites_easy));
		
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

		writeFile("./myOut193_Wh.txt",wmnohfm1+" "+wmnolfm1+" "+wmnohfm2+" "+wmnolfm2,append = T);
		writeFile("./myOut193_Oc.txt",omnohfm1+" "+omnolfm1+" "+omnohfm2+" "+omnolfm2,append = T);
		writeFile("./myOut193_Fr.txt",fmnohfm1+" "+fmnolfm1+" "+fmnohfm2+" "+fmnolfm2,append = T);
			
			
	}
	if(sim.generation == END)
	{
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
