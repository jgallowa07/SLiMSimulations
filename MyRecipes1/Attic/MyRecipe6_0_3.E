//This model has mutations for both freshwater and marine that are recessive (complete) and dominant (complete)

initialize() {

	pathToOutputFile = "/projects/haldane/shared/SLiMSimulations/Output1/MyRecipe6_0_3/";
	//pathToOutputFile = "~/Documents/ResearchStuff/CODEBASS/Assigned/Output1/MyRecipe6_0_3Local/";

	//KNOBS//

	endOfSimulation = 40000;

	samplingInterval = 100;

	numberOfLakes = 10;	
	
	optimumPhenoMarine = 10;
	
	optimumPhenoFresh = -10;

	individualsInOcean = 2000;
	
	individualsInFresh = 2000;

	numberOfEffectRegions = 10;

	sizeOfEffectRegions = 100;

	sizeOfNonEffectBuffer = 4950;
	
	mutationRate = 1e-7;
		
	recombinationRate = 1e-6;

	introduce = 20000;

	migrationRate = 0.005;
		

	individualsPerFresh = individualsInFresh/numberOfLakes;	
	lengthOfPops = asFloat(numberOfLakes);
	numberOfSamplingPoints = endOfSimulation/samplingInterval;	
	
	defineConstant("END",endOfSimulation);
	defineConstant("INT",samplingInterval);
	defineConstant("NLAKES",numberOfLakes);
	defineConstant("IIO", individualsInOcean);
	defineConstant("IIF", individualsInFresh);	
	defineConstant("LOP", lengthOfPops);
	defineConstant("IPF",individualsPerFresh);
	defineConstant("PATH",pathToOutputFile);
	defineConstant("MR", mutationRate);
	defineConstant("C", numberOfEffectRegions);    	// number of Effect Regions;
	defineConstant("W", sizeOfNonEffectBuffer);  	// size of neutral buffer on each side
	defineConstant("P", sizeOfEffectRegions);   	// Size of contribution (Effect) Loci 
	defineConstant("INP", introduce);
	defineConstant("MGR", migrationRate);

	initializeSLiMOptions(dimensionality="x");
	initializeMutationRate(MR); 	

	initializeMutationType("m1", 0.5, "g", 0.5, 1.0);  //Marine Dominant 
	initializeMutationType("m2", 0.5, "g", -0.5, 1.0); //Freshwater Dominant
	initializeMutationType("m3", 0.5, "f", 0.0);       //Nuetral Benefitial Mutations
	initializeMutationType("m4", 0.5, "g", 0.5, 1.0);  //Marine Recessive
	initializeMutationType("m5", 0.5, "g", -0.5, 1.0); //Freshwater Recessive
	initializeMutationType("m6", 0.5, "g", 0.6, 1.0);  //Marine Additive
	initializeMutationType("m7", 0.5, "g", -0.6, 1.0); //Freshwater Additive

	initializeGenomicElementType("g1", m3, 1.0);                       	//NonEffect Genomic Element
	initializeGenomicElementType("g2", c(m1,m2,m3,m4,m5,m6,m7), c(1.0,1.0,1.0,1.0,1.0,1.0,1.0)); 	//Effect Genomic Element

	m1.mutationStackPolicy = "l";
	m2.mutationStackPolicy = "l";
	m3.mutationStackPolicy = "l";
	m4.mutationStackPolicy = "l";
	m5.mutationStackPolicy = "l";
	m6.mutationStackPolicy = "l";
	m7.mutationStackPolicy = "l";
	
	//initialize genomic arch.
	// |----------NonEffect-----------|-Effect(g3)-|----------NonEffect----------|

	pos = 0;
	q = NULL;
	for (i in 1:C)
	{
		initializeGenomicElement(g1, pos, pos + W-1);    
		pos = pos + W;

		initializeGenomicElement(g2, pos, pos + P-1);
		//section = pos:(pos+P-1);
		//q = c(q,section);
		q = c(q,c(pos,pos+P-1));
		pos = pos + P;	

		initializeGenomicElement(g1, pos, pos + W-1);
		pos = pos + W;
	}

	catn("end of chromosome position: "+pos);
	totalNumOfLoci = pos;

	defineConstant("TL", totalNumOfLoci);
	defineConstant("Q", q);     // remember our QTL positions
		
	catn(Q);	
	catn(TL);

	initializeRecombinationRate(recombinationRate);	

	//Interation Function for Freshwater mateChoice() callback...

	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.02);
	i1.setInteractionFunction("n", 0.5, 0.02);

	//Interaction Function for Oceanic mateChoice() callback...

	initializeInteractionType(2, "x", reciprocal=T, maxDistance=1.0);
	i2.setInteractionFunction("n", 0.5, 0.5);
	
	writeFile(PATH+"AveragePheno.txt",numberOfSamplingPoints+" "+INT+" "+((INP/INT)-1),append = F);
	writeFile(PATH+"MeanVarience.txt",numberOfSamplingPoints+" "+INT,append = F);
	//writeFile(PATH+"AveragePheno.txt",numberOfSamplingPoints+" "+INT+" "+NLAKES,append = F);
	//writeFile(PATH+"myData.ped"," ",append = F);
	//writeFile(PATH+"AverageFst.txt",numberOfSamplingPoints+" "+INT,append = F);
	//writeFile(PATH+"MatrixAllelesAll.txt",100+" "+TL,append = F);
	//writeFile(PATH+"PhenoDist.txt",100+" "+TL,append = F);
	//writeFile(PATH+"Correlations.txt",asString(TL),append = F);
	
}


1 late(){
	sim.addSubpop("p1", IIO);  //Marine 
	sim.addSubpop("p2", IIF);  //Freshwater
	
	//set length of pops equal to 1 unit per lake

	p1.setSpatialBounds(c(0.0, LOP));
	p2.setSpatialBounds(c(0.0, LOP));

	//set migration to 1 individual per generation from each pop to another 

	//migrationRateIntoOcean = asFloat(NLAKES/IIF);
	//migrationRateIntoFresh = asFloat(1/IIO);
	
	p1.setMigrationRates(p2,MGR);
	p2.setMigrationRates(p1,MGR);
	
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


20000 late(){

	sim.addSubpopSplit("p3",IIF,p1);
	p1.setSpatialBounds(c(0.0, LOP));
	
	lakeIndividuals2 = p3.individuals;
	
	migrationRateIntoOcean = asFloat(NLAKES/IIF)/2.0;
	//migrationRateIntoFresh = asFloat(1/IIO);
	
	p1.setMigrationRates(p3,MGR);
	p3.setMigrationRates(p1,MGR);

	for (i in 0:(NLAKES-1)){
		for (j in 0:(IPF-1)){
			lakeIndividuals2[IPF*i + j].setSpatialPosition((i+1)-(0.5));						
		}
	}
}


1:19999 late() {  //ALMOST ALL SAMPLING
	i1.evaluate();
	i2.evaluate();		
	Oceanic_mean_phenotype = NULL;
	Freshwater_mean_phenotype = NULL;
	substitutionOffset  = 0;

	for (individual in sim.subpopulations.individuals) //Calculate fitness (.tagF) 
	{
		dominantMarine = 2 * sum(individual.uniqueMutationsOfType(m1).selectionCoeff);
		dominantFreshwater = 2 * sum(individual.uniqueMutationsOfType(m2).selectionCoeff);

		aM = sum(individual.genomes.mutationsOfType(m4).selectionCoeff);
		aF = sum(individual.genomes.mutationsOfType(m5).selectionCoeff);
		dM = 2 * sum(individual.uniqueMutationsOfType(m4).selectionCoeff);
		dF = 2 * sum(individual.uniqueMutationsOfType(m5).selectionCoeff);

		recessiveMarine = (2 * aM) - dM; 
		recessiveFreshwater = (2 * aF) - dF;

		additiveMarine = sum(individual.genomes.mutationsOfType(m6).selectionCoeff);
		additiveFreshwater = sum(individual.genomes.mutationsOfType(m7).selectionCoeff);
	
		individual.tagF = dominantMarine + recessiveFreshwater + dominantFreshwater + recessiveMarine + additiveMarine + additiveFreshwater;		
	}

	
	if(sim.generation == 19999){

	//get the standard deviation HERE	
	
		phenotypes = rep(0.0,IIO);
		count = 0;
		for(individual in p1.individuals){
			phenotypes[count] = individual.tagF;
			count = count + 1;
		}
		standard_deviation = sd(phenotypes);
	
		writeFile("/projects/haldane/shared/SLiMSimulations/Output1/SD.txt","6_0_3: SD with MiGRate "+MGR+" is: "+standard_deviation,append = T);
		//writeFile("~/Documents/ResearchStuff/CODEBASS/Assigned/Output1/SD.txt","SD with MiGRate "+MGR+" is: "+standard_deviation,append = T);
		//catn(standard_deviation);
			
	}	

	//The entire if{} block below is for sampling data at a given interval, INT	
	
	if((sim.generation % INT == 0) | (sim.generation == END))
	{
		//TODO Number of mutations per individual. 
	
		cat("-------------------------------\n");
		cat("Output1 for generation " + sim.generation + ":\n\n");
		Oceanic_mean_phenotype =  mean(p1.individuals.tagF);
		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
		catn("Oceanic mean pheno: "+Oceanic_mean_phenotype);
		catn("Freshwater mean pheno: "+Freshwater_mean_phenotype);			

		writeFile(PATH+"AveragePheno.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype,append = T);

		allMuts = sim.mutations;
	
		numCells = 1000;
		cellSize = TL/numCells;
		lociPositions = 0:(TL-1);
		cellCounts = rep(0,numCells);
		breakPoints = lociPositions[lociPositions % cellSize == 0];
		index = 0;
		totalNumIndividuals = (IIO) + (2*IIF);

		for (i in breakPoints){

			mutsInCell = allMuts[(allMuts.position >= i) & (allMuts.position < (i+cellSize))];
			cellCounts[index] = size(mutsInCell);

		//  --number of muts across population--
		//	mutFreqs = sim.mutationFrequencies(NULL,mutsInCell);
		//	totalNumMuts = sum(mutFreqs * totalNumIndividuals);
		//	cellCounts[index] = asInteger(totalNumMuts);

			index = index + 1;

		}

		catn(asString(cellCounts));				

		Mean = mean(cellCounts);
		catn("Mean: " + asString(Mean));
	
		z = sd(cellCounts);	
		Varience = z*z;
		catn("Varience: " + asString(Varience));

		writeFile(PATH+"MeanVarience.txt",Mean+" "+Varience,append = T);

	}

	//Set tag equal to a value to be returned by fitness callback 
	
	for (individual in p2.individuals){
		individual.tagF = dnorm(-10.0 - individual.tagF,0.0,15.0);
	}
	for (individual in p1.individuals){
		individual.tagF = dnorm(10 - individual.tagF,0.0,15.0);
	}

	//Normalize the fitness in the lakes so an equivelent mean number of parents get chosen from each.
 
	for (l in 1:NLAKES){
		locationOfLake = (l - 0.5);
		lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
		sumOfPhenos = sum(lake_individuals.tagF);
		lake_individuals.tagF =(lake_individuals.tagF*IPF)/(sumOfPhenos);
	}
	
}


20000: late(){

	i1.evaluate();
	i2.evaluate();		
	Oceanic_mean_phenotype = NULL;
	Freshwater_mean_phenotype = NULL;
	substitutionOffset  = 0;

	for (individual in sim.subpopulations.individuals) //Calculate fitness (.tagF) 
	{
		dominantMarine = 2 * sum(individual.uniqueMutationsOfType(m1).selectionCoeff);
		dominantFreshwater = 2 * sum(individual.uniqueMutationsOfType(m2).selectionCoeff);

		aM = sum(individual.genomes.mutationsOfType(m4).selectionCoeff);
		aF = sum(individual.genomes.mutationsOfType(m5).selectionCoeff);
		dM = 2 * sum(individual.uniqueMutationsOfType(m4).selectionCoeff);
		dF = 2 * sum(individual.uniqueMutationsOfType(m5).selectionCoeff);

		recessiveMarine = (2 * aM) - dM; 
		recessiveFreshwater = (2 * aF) - dF;

		additiveMarine = sum(individual.genomes.mutationsOfType(m6).selectionCoeff);
		additiveFreshwater = sum(individual.genomes.mutationsOfType(m7).selectionCoeff);
	
		individual.tagF = dominantMarine + recessiveFreshwater + dominantFreshwater + recessiveMarine + additiveMarine + additiveFreshwater;		
	
	}

	
	if((sim.generation % INT == 0) | (sim.generation == END))
	{
		//TODO Number of mutations per individual. 
	
		//---------------------Sample Average Phenotype----------------------------------------
		cat("-------------------------------\n");
		cat("Output1 for generation " + sim.generation + ":\n\n");
		Oceanic_mean_phenotype =  mean(p1.individuals.tagF);
		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
		newFreshwater_mean_phenotype = mean(p3.individuals.tagF);
		catn("Oceanic mean pheno: "+Oceanic_mean_phenotype);
		catn("Freshwater mean pheno: "+Freshwater_mean_phenotype);	
		catn("new Freshwater mean pheno: "+newFreshwater_mean_phenotype);		

		writeFile(PATH+"AveragePheno.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype+" "+newFreshwater_mean_phenotype,append = T);

		allMuts = sim.mutations;
	
		numCells = 1000;
		cellSize = TL/numCells;
		lociPositions = 0:(TL-1);
		cellCounts = rep(0,numCells);
		breakPoints = lociPositions[lociPositions % cellSize == 0];
		index = 0;
		totalNumIndividuals = (IIO) + (2*IIF);

		for (i in breakPoints){

			mutsInCell = allMuts[(allMuts.position >= i) & (allMuts.position < (i+cellSize))];
			cellCounts[index] = size(mutsInCell);

		//  --number of muts across population--
		//	mutFreqs = sim.mutationFrequencies(NULL,mutsInCell);
		//	totalNumMuts = sum(mutFreqs * totalNumIndividuals);
		//	cellCounts[index] = asInteger(totalNumMuts);

			index = index + 1;

		}

		catn(asString(cellCounts));				

		Mean = mean(cellCounts);
		catn("Mean: " + asString(Mean));
	
		z = sd(cellCounts);	
		Varience = z*z;
		catn("Varience: " + asString(Varience));

		writeFile(PATH+"MeanVarience.txt",Mean+" "+Varience,append = T);
		


	}

	if(sim.generation == END)
	{

//		m1s = sim.mutationsOfType(m1);			
//		m2s = sim.mutationsOfType(m2);			
//		m4s = sim.mutationsOfType(m4);			
//		m5s = sim.mutationsOfType(m5);			
//		m6s = sim.mutationsOfType(m6);			
//		m7s = sim.mutationsOfType(m7);			

//		allEffectMuts = c(m1s,m2s,m4s,m5s,m6s,m7s);



		sim.simulationFinished();
	}
	for (individual in p3.individuals){
		individual.tagF = dnorm(-10.0 - individual.tagF,0.0,15.0);
	}
	for (individual in p2.individuals){
		individual.tagF = dnorm(-10.0 - individual.tagF,0.0,15.0);
	}
	for (individual in p1.individuals){
		individual.tagF = dnorm(10 - individual.tagF,0.0,15.0);
	}

	//Normalize the fitness in the lakes so an equivelent mean number of parents get chosen from each.
 
	for (l in 1:NLAKES){
		locationOfLake = (l - 0.5);
		lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
		sumOfPhenos = sum(lake_individuals.tagF);
		lake_individuals.tagF =(lake_individuals.tagF*IPF)/(sumOfPhenos);
		lake_individuals2 = p3.individuals[abs(p3.individuals.spatialPosition - locationOfLake) < 1e-8];
		sumOfPhenos2 = sum(lake_individuals2.tagF);
		lake_individuals2.tagF =(lake_individuals2.tagF*IPF)/(sumOfPhenos2);
	}

}

fitness(m1) {return 1.0;}
fitness(m2) {return 1.0;}
fitness(m4) {return 1.0;}
fitness(m5) {return 1.0;}
fitness(m6) {return 1.0;}
fitness(m7) {return 1.0;}
fitness(NULL, p1) {
	return individual.tagF;
}
fitness(NULL, p2) {
	return individual.tagF;
}
20000: fitness(NULL, p3){
	return individual.tagF;
}
1: mateChoice(p2) {
	return  i1.strength(individual) * weights; 
}
1: mateChoice(p1){
	return i2.strength(individual) * weights;
}
20000: mateChoice(p3){
	return i1.strength(individual) * weights;
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
20000: modifyChild(p3){

	if(sourceSubpop != subpop){  //migrant from Lake -> Oceanic 
		pos = parent1.spatialPosition + rnorm(1,0,0.5);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}else{
		child.setSpatialPosition(parent1.spatialPosition);
	}
	return T;
}



40001 late() {sim.simulationFinished();}
