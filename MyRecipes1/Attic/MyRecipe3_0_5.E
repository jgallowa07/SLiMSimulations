initialize() {



	pathToOutputFile = "/projects/haldane/SLiMSimulations/Output1/MyRecipe3_0_5/";
	pathToOutputFile = "~/Documents/ResearchStuff/CODEBASS/Assigned/Output1/MyRecipe3_0_5Local/";

	//KNOBS//

	endOfSimulation = 100;

	samplingInterval = 20;

	numberOfLakes = 10;	

	freshwaterOptimumFitness = -10;
	
	marineOptimumFitness = 10;

	individualsInOcean = 2000;
	
	individualsInFresh = 2000;

	numberOfEffectRegions = 10;

	sizeOfEffectRegions = 100;

	sizeOfNonEffectBuffer = 4950;
	
	mutationRate = 1e-7;
		
	recombinationRate = 1e-3;
	

	individualsPerFresh = individualsInFresh/numberOfLakes;	
	lengthOfPops = asFloat(numberOfLakes);
	numberOfSamplingPoints = endOfSimulation/samplingInterval + 1;	
	
	defineConstant("END",endOfSimulation);
	defineConstant("INT",samplingInterval);
	defineConstant("NLAKES",numberOfLakes);
	defineConstant("FO", freshwaterOptimumFitness);
	defineConstant("MO", marineOptimumFitness);
	defineConstant("IIO", individualsInOcean);
	defineConstant("IIF", individualsInFresh);	
	defineConstant("LOP", lengthOfPops);
	defineConstant("IPF",individualsPerFresh);
	defineConstant("PATH",pathToOutputFile);
	defineConstant("MR", mutationRate);
	defineConstant("C", numberOfEffectRegions);    	// number of Effect Regions;
	defineConstant("W", sizeOfNonEffectBuffer);  	// size of neutral buffer on each side
	defineConstant("P", sizeOfEffectRegions);   	// Size of contribution (Effect) Loci 
	

	initializeSLiMOptions(dimensionality="x");
	initializeMutationRate(MR); 	

	initializeMutationType("m1", 0.5, "g", 0.5, 1.0);  //Oceanic Benefitial Mutations
	initializeMutationType("m2", 0.5, "g", -0.5, 1.0); //Freshwater Benefictial Mutations
	initializeMutationType("m3", 0.5, "f", 0.0);       //Nuetral Benefitial Mutations
	initializeMutationType("m4", 0.5, "f", 0.0);       //deleterious Mutations
	
	initializeGenomicElementType("g1", m3, 1.0);                       	//NonEffect Genomic Element
	initializeGenomicElementType("g2", c(m1,m2,m3), c(1.0,1.0,1.0)); 	//Effect Genomic Element

	m1.mutationStackPolicy = "l";
	m2.mutationStackPolicy = "l";
	
	//initialize genomic arch.
	// |----------NonEffect-----------|-Effect-|----------NonEffect----------|	  Repeat for number of Effect Regions ...

	pos = 0;
	q = NULL;
	
	for (i in 1:C)
	{
		initializeGenomicElement(g1, pos, pos + W-1);    
		pos = pos + W;
	
		initializeGenomicElement(g2, pos, pos + P-1);
		section = pos:(pos+P-1);
		q = c(q,section);
		pos = pos + P;	

		initializeGenomicElement(g1, pos, pos + W-1);
		pos = pos + W;
	}

	catn("end of chromosome position: "+pos);
	totalNumOfLoci = pos;

	defineConstant("TL", totalNumOfLoci);
	defineConstant("Q", q);     // remember our QTL positions

	initializeRecombinationRate(recombinationRate);	

	//Interation Function for Freshwater mateChoice() callback...

	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.02);
	i1.setInteractionFunction("n", 0.5, 0.02);

	//Interaction Function for Oceanic mateChoice() callback...

	initializeInteractionType(2, "x", reciprocal=T, maxDistance=1.0);
	i2.setInteractionFunction("n", 0.5, 0.5);
	
	writeFile(PATH+"AveragePheno.txt",numberOfSamplingPoints+" "+INT+" "+NLAKES,append = F);
	writeFile(PATH+"AverageFst.txt",numberOfSamplingPoints+" "+INT,append = F);
	writeFile(PATH+"MatrixAllelesAll.txt",100+" "+TL,append = F);
	writeFile(PATH+"PhenoDist.txt",100+" "+TL,append = F);
	
}


1 late(){
	sim.addSubpop("p1", IIO);  //Marine 
	sim.addSubpop("p2", IIF);  //Freshwater
	
	//set length of pops equal to 1 unit per lake

	p1.setSpatialBounds(c(0.0, LOP));
	p2.setSpatialBounds(c(0.0, LOP));

	//set migration to 1 individual per generation from each pop to another 

	migrationRate = asFloat(NLAKES/IIF);
	
	p1.setMigrationRates(p2,migrationRate);
	p2.setMigrationRates(p1,migrationRate);
	
	//unifomily distribute oceanic individuals.

	for (ind in p1.individuals){
		ind.setSpatialPosition(p1.pointUniform());
	}
	
	lakeIndividuals = p2.individuals;

	//Put an even amount of individuals in each "lake" (location)

	for (i in 0:(NLAKES-1)){
		for (j in 0:(IPF-1)){
			lakeIndividuals[IPF*i + j].setSpatialPosition((i+1)-(0.5));						
		}
	}

}


1: late() {  //ALMOST ALL SAMPLING
	i1.evaluate();
	i2.evaluate();		
	Oceanic_mean_phenotype = NULL;
	Freshwater_mean_phenotype = NULL;
	substitutionOffset  = 0;

	for (fix in sim.substitutions){
		substitutionOffset = substitutionOffset + fix.selectionCoeff;
	}	

	for (individual in sim.subpopulations.individuals) //Calculate fitness (.tagF) 
	{
		EffectMuts = c(individual.genomes.mutationsOfType(m1),individual.genomes.mutationsOfType(m2));
		individual.tagF = substitutionOffset + sum(EffectMuts.selectionCoeff);		
	}
	
	if((sim.generation == END/2) | (sim.generation == END))
	{
			
		//TODO: PUT MATRIX CALC HERE FOR MULIPLE SAMPLING TIMES THROUGH SIMUALTION

	}
	

	//The entire if{} block below is for sampling data at a given interval, INT	
	
	if((sim.generation  % INT == 0) | (sim.generation == END))
	{
		//TODO Number of mutations per individual. 
	
		cat("-------------------------------\n");
		cat("Output1 for generation " + sim.generation + ":\n\n");
		Oceanic_mean_phenotype =  mean(p1.individuals.tagF);
		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
		catn("Oceanic mean pheno: "+Oceanic_mean_phenotype);
		catn("Freshwater mean pheno: "+Freshwater_mean_phenotype);			


		//---------------------Sample Average Phenotype----------------------------------------
		
		lake_phenos = NULL;		

		for (i in 1:NLAKES)
		{
			//Pheno sampling	
			
			locationOfLake = (i - 0.5);
			lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
			lake_mean_pheno = mean(lake_individuals.tagF);
			lake_phenos = c(lake_phenos,lake_mean_pheno);			
		}
		
	
		lake_phenos_string = "";
		
		for (i in 0:(NLAKES-1))
		{
			lake_phenos_string = lake_phenos_string+asString(lake_phenos[i])+" ";
		}
			
		writeFile(PATH+"AveragePheno.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype+" "+lake_phenos_string,append = T);

		//-----------------------Between Ocean/Freshwater Fst Sample-----------------------------

		p1_p = sim.mutationFrequencies(p1);
		p2_p = sim.mutationFrequencies(p2);
		mean_p = (p1_p + p2_p) / 2.0;
		pq = (p1_p * (1 - p1_p) + p2_p * (1 - p2_p))/2.0;
		mean_pq = mean_p * (1 - mean_p);
		fst = 1 - (pq / mean_pq);		
		meanOceanFreshFst = mean(fst);
		//catn("mean Ocean/Fresh Fst: "+ meanOceanFreshFst);
		
		//-----------------------Between Ocean/Freshwater Fst Sample (Neutral Muts)-----------------------------

		neutralMutations = sim.mutationsOfType(m3);

		p1_n = sim.mutationFrequencies(p1,neutralMutations);
		p2_n = sim.mutationFrequencies(p2,neutralMutations);
		mean_n = (p1_n + p2_n) / 2.0;
		pqn = (p1_n * (1 - p1_n) + p2_n * (1 - p2_n))/2.0;
		mean_pqn = mean_n * (1 - mean_n);
		fst_n = 1 - (pqn / mean_pqn);	
		
		meanOceanFreshFstNeutral = mean(fst_n);
		//catn("mean Ocean/Fresh Fst (Neutral mutations only): "+ meanOceanFreshFstNeutral);
		
		//-----------------------Between Ocean/Freshwater Fst Sample (Effect Muts)-----------------------------
		
		m1mutations = sim.mutationsOfType(m1);
		m2mutations = sim.mutationsOfType(m2);
	
		effectMutations = c(m1mutations,m2mutations);

		p1_e = sim.mutationFrequencies(p1,effectMutations);
		p2_e = sim.mutationFrequencies(p2,effectMutations);
		mean_e = (p1_e + p2_e) / 2.0;
		pqe = (p1_e * (1 - p1_e) + p2_e * (1 - p2_e))/2.0;
		mean_pqe = mean_e * (1 - mean_e);
		fst_e = 1 - (pqe / mean_pqe);	
		
		meanOceanFreshFstEffect = mean(fst_e);
		//catn("mean Ocean/Fresh Fst (Effect mutations only): "+ meanOceanFreshFstEffect);
		
		
			 
		//-----------------------Between Lakes Fst Sample-----------------------------
		
		muts = sim.mutations;
		freqs = rep(0,size(muts));
		mean_p = freqs;
		mean_het = freqs;
	
		for (i in 1:NLAKES)
		{	
			locationOfLake = (i - 0.5);
			lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
			
			freqs = rep(0, size(muts));
			
			for (geno in lake_individuals.genomes)
			{
				freqs = freqs + asInteger(geno.containsMutations(muts));		
			}
			freqs = freqs / (2*size(lake_individuals));
			mean_p = mean_p + freqs;
			mean_het = mean_het + (freqs * (1-freqs));
		} 
		
		mean_p = mean_p/NLAKES;        
		temp = mean_p;
		mean_het = mean_het/NLAKES;         

		mean_p = mean_p[mean_p != 0 & mean_het != 0];
		mean_het = mean_het[temp != 0 & mean_het != 0];		
	
		between_fst = mean(1 - mean_het / (mean_p * (1 - mean_p)));
				
		//catn("between lakes Fst: "+ between_fst);							


		//-----------------------Between Lakes Fst Sample (Neutral Mutations)-----------------------------

		
		freqs_n = rep(0,size(neutralMutations));
		mean_p_n = freqs_n;
		mean_het_n = freqs_n;
	
		for (i in 1:NLAKES)
		{	
			locationOfLake = (i - 0.5);
			lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
			
			freqs_n = rep(0, size(neutralMutations));
			
			for (geno in lake_individuals.genomes)
			{
				freqs_n = freqs_n + asInteger(geno.containsMutations(neutralMutations));		
			}
			freqs_n = freqs_n / (2*size(lake_individuals));
			mean_p_n = mean_p_n + freqs_n;
			mean_het_n = mean_het_n + (freqs_n * (1-freqs_n));
		} 
		
		mean_p_n = mean_p_n/NLAKES;        
		temp_n = mean_p_n;
		mean_het_n = mean_het_n/NLAKES;         

		mean_p_n = mean_p_n[mean_p_n != 0 & mean_het_n != 0];
		mean_het_n = mean_het_n[temp_n != 0 & mean_het_n != 0];		
	
		between_fst_n = mean(1 - mean_het_n / (mean_p_n * (1 - mean_p_n)));
				
		//catn("between lakes Fst (neutral): "+ between_fst_n);							

		//-----------------------Between Lakes Fst Sample (Effect Mutations)-----------------------------
		
		freqs_e = rep(0,size(effectMutations));
		mean_p_e = freqs_e;
		mean_het_e = freqs_e;
	
		for (i in 1:NLAKES)
		{	
			locationOfLake = (i - 0.5);
			lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
			
			freqs_e = rep(0, size(effectMutations));
			
			for (geno in lake_individuals.genomes)
			{
				freqs_e = freqs_e + asInteger(geno.containsMutations(effectMutations));		
			}
			freqs_e = freqs_e / (2*size(lake_individuals));
			mean_p_e = mean_p_e + freqs_e;
			mean_het_e = mean_het_e + (freqs_e * (1-freqs_e));
		} 
		
		mean_p_e = mean_p_e/NLAKES;        
		temp_e = mean_p_e;
		mean_het_e = mean_het_e/NLAKES;         

		//catn("mean_p_e: " + mean_p_e);
		//catn("mean_het_e: " + mean_het_e);

		mean_p_e = mean_p_e[mean_p_e != 0 & mean_het_e != 0];
		mean_het_e = mean_het_e[temp_e != 0 & mean_het_e != 0];		
	
		//catn(mean_p_e);
		//catn(mean_het_e);

		between_fst_e = mean(1 - mean_het_e / (mean_p_e * (1 - mean_p_e)));
				
		//catn("between lakes Fst (effect): "+ between_fst_e);							

	
		//-----------------------Between Ocean Fst Sample-----------------------------

		//muts = sim.mutations;
		freqs = rep(0,size(muts));
		mean_p = freqs;
		mean_het = freqs;
	
		for (i in 1:NLAKES)
		{	
			leftBound = (i - 1);
			rightBound = i;

			region_individuals = p1.individuals[p1.individuals.x >= leftBound & p1.individuals.x < rightBound];
						
			freqs = rep(0, size(muts));
			
			for (geno in region_individuals.genomes)
			{
				freqs = freqs + asInteger(geno.containsMutations(muts));		
			}
			freqs = freqs / (2*size(region_individuals));
			mean_p = mean_p + freqs;
			mean_het = mean_het + (freqs * (1-freqs));
		} 
		
		mean_p = mean_p/NLAKES;        
		temp = mean_p;
		mean_het = mean_het/NLAKES;         

		mean_p = mean_p[mean_p != 0 & mean_het != 0];
		mean_het = mean_het[temp != 0 & mean_het != 0];	
	
		between_marine_fst = mean(1 - mean_het / (mean_p * (1 - mean_p)));
				
		//catn("between Marine Regions Fst: "+ between_marine_fst);						
			
		//writeFile(PATH+"AverageFst.txt",meanOceanFreshFst+" "+between_fst+" "+between_marine_fst+" "+meanOceanFreshFstNeutral+" "+meanOceanFreshFstEffect+" "+between_fst_n+" "+between_fst_e,append = T);
	
	}
	if(sim.generation == END)
	{
	
		//------------------------------------------Mutations Effect, Frequency, Position ------------------------------
		
		muts = sim.mutations;
		p1muts = unique(p1.genomes.mutations);
		p2muts = unique(p2.genomes.mutations);

		num = size(muts);
		p1num = size(p1muts);
		p2num = size(p2muts);

		catn("whole: "+num);
		catn("Oceanic: "+p1num);
		catn("Freshwater: "+p2num);	

		mutfreq = sim.mutationFrequencies(NULL,muts);
		muteff = muts.selectionCoeff;
		mutpos = muts.position;
		writeFile(PATH+"WholeMutations.txt",num+" ",append = F);
		for (i in 0:(num-1))
		{
			writeFile(PATH+"WholeMutations.txt",mutfreq[i]+" "+mutpos[i]+" "+muteff[i],append = T);
		}

		mutfreq1 = sim.mutationFrequencies(p1,p1muts);
		muteff1 = p1muts.selectionCoeff;
		mutpos1 = p1muts.position;
		writeFile(PATH+"OceanicMutations.txt",p1num+" ",append = F);
		for (i in 0:(p1num-1))
		{
			writeFile(PATH+"OceanicMutations.txt",mutfreq1[i]+" "+mutpos1[i]+" "+muteff1[i],append = T);
		}

		mutfreq2 = sim.mutationFrequencies(p2,p2muts);
		muteff2 = p2muts.selectionCoeff;
		mutpos2 = p2muts.position;
		writeFile(PATH+"FreshwaterMutations.txt",p2num+" ",append = F);
		for (i in 0:(p2num-1))
		{
			writeFile(PATH+"FreshwaterMutations.txt",mutfreq2[i]+" "+mutpos2[i]+" "+muteff2[i],append = T);
		}
		
		//sample a single individual haplotype from each population.
		
		oceanic_geno = sample(p1.individuals.genomes,1);
		lake_geno = sample(p2.individuals.genomes,1);
		
		oceanic_ind_muts = c(oceanic_geno.mutationsOfType(m1),oceanic_geno.mutationsOfType(m2));
		num_Oind_muts = size(oceanic_ind_muts);
			
		oceanic_ind_mut_freq = sim.mutationFrequencies(p1,oceanic_ind_muts);
		oceanic_ind_mut_eff = oceanic_ind_muts.selectionCoeff;
		oceanic_ind_mut_pos = oceanic_ind_muts.position;
	

		lake_ind_muts = c(lake_geno.mutationsOfType(m1),lake_geno.mutationsOfType(m2));
		num_Lind_muts = size(lake_ind_muts);
					
		lake_ind_mut_freq = sim.mutationFrequencies(p2,lake_ind_muts);
		lake_ind_mut_eff = lake_ind_muts.selectionCoeff;
		lake_ind_mut_pos = lake_ind_muts.position;

		writeFile(PATH+"IndOceanicMutations.txt",num_Oind_muts+" ",append = F);
		for (i in 0:(num_Oind_muts-1))
		{
			//writeFile(PATH+"IndOceanicMutations.txt",oceanic_ind_mut_freq[i]+" "+oceanic_ind_mut_pos[i]+" "+oceanic_ind_mut_eff[i],append = T);
		}
		
		writeFile(PATH+"IndFreshwaterMutations.txt",num_Lind_muts+" ",append = F);
		for (i in 0:(num_Lind_muts-1))
		{
			//writeFile(PATH+"IndFreshwaterMutations.txt",lake_ind_mut_freq[i]+" "+lake_ind_mut_pos[i]+" "+lake_ind_mut_eff[i],append = T);
		}
		//-----------------------------------------Fst by Mutation------------------------------------------------


		//Sample Fst between Ocean and Freshwater combined.
		
		p1_p = sim.mutationFrequencies(p1,muts); //get frequencies of all mutations in terms of p1
		p2_p = sim.mutationFrequencies(p2,muts); //get frequencies of all mutations in terms of p2
		mean_p = (p1_p + p2_p) / 2.0;       //average the frequencies. 
		pq = (p1_p * (1 - p1_p) + p2_p * (1 - p2_p))/2.0;
		mean_pq = mean_p * (1 - mean_p);
		fst = 1 - (pq / mean_pq);

		catn("fstsize: "+ size(fst));
		
		writeFile(PATH+"OceanFreshwaterFst.txt",asString(size(fst)),append = F);
		for (i in 0:(num-1))
		{
			writeFile(PATH+"OceanFreshwaterFst.txt",fst[i]+" "+mutpos[i] ,append = T);
		}
		
	
		//-----------------------------------------Matrix-Alleles-------------------------------------------------

		effectPos = effectMutations.position;
		neutralPos = neutralMutations.position;		
		
	
		//catn("sizes: "+ size(muts) + " "+ size(mutPositions) + " " + size(marineCommon) );
	
		sampleGenomes = NULL;
	
		for (i in 1:NLAKES)
		{	
			leftBound = (i - 1);
			rightBound = i;

			region_individuals = p1.individuals[p1.individuals.x >= leftBound & p1.individuals.x < rightBound];

			sam = sample(region_individuals.genomes,10);
	
			sampleGenomes = c(sampleGenomes,sam);
		}	
	
		for (i in 1:NLAKES)
		{		
			locationOfLake = (i - 0.5);
			lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
		
			sam = sample(lake_individuals.genomes,10);
			sampleGenomes = c(sampleGenomes,sam);				
		}
		
		count = 1;		
	
		for (samp in sampleGenomes)
		{
			containsEffect = samp.containsMutations(effectMutations);
			containsNeutral = samp.containsMutations(neutralMutations);
	
			sampEffPos = effectPos[containsEffect];
			sampNeuPos = neutralPos[containsNeutral];	
		
			for (a in sampEffPos)
			{
				writeFile(PATH+"MatrixAllelesAll.txt",count + " " + asString(a) + " " + asString(-1), append = T);
			}
			
			for (b in sampNeuPos)
			{
				writeFile(PATH+"MatrixAllelesAll.txt",count + " " + asString(b) + " " + asString(1), append = T);
			}
	
		
			count = count + 1;	
	
		}	
	
	
		//----------------------------------- Phenotype Histogram -------------------------------------------

		for(ind in sim.subpopulations.individuals)
		{
			writeFile(PATH+"PhenoDist.txt",asString(ind.tagF),append = T);
		}
		
		//----------------------------------- Simulation Finished -------------------------------------------

		sim.simulationFinished();

	}

	//Set tag equal to a value to be returned by fitness callback 
	
	for (individual in p2.individuals){
		individual.tagF = dnorm(FO - individual.tagF,0.0,15.0);
	}
	for (individual in p1.individuals){
		individual.tagF = dnorm(MO - individual.tagF,0.0,15.0);
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
	//numDeleterious = individual.countOfMutationsOfType(m4);
	//effect = numDeleterious/10.0;
	return individual.tagF;
}
fitness(NULL, p2) {
	//numDeleterious = individual.countOfMutationsOfType(m4);
	//effect = numDeleterious/10.0;
	return individual.tagF;
}
1: mateChoice(p2) {
	return  i1.strength(individual) * weights; 
}
1: mateChoice(p1){
	return i2.strength(individual) * weights;
}

modifyChild(p1){

	if(sourceSubpop != subpop){  //mirant from the Ocean -> Lake 
		num = parent1.spatialPosition + rnorm(1,0,0.5);
		num = p2.pointReflected(num);
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

	if(sourceSubpop != subpop){  //migrant from Lake -> Oceanic 
		pos = parent1.spatialPosition + rnorm(1,0,0.5);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}else{
		child.setSpatialPosition(parent1.spatialPosition);
	}
	return T;
}


160000 late() {sim.simulationFinished();}
