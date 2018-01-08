//This model has mutations for both freshwater and marine that are recessive (complete) and dominant (complete)

initialize() {

	pathToOutputFile = "/projects/haldane/shared/SLiMSimulations/Output1/MyRecipe8_0_5/";
	//pathToOutputFile = "~/Documents/ResearchStuff/CODEBASS/Assigned/Output1/MyRecipe8_0_5Local/";

	//KNOBS//

	endOfSimulation = 80000;

	samplingInterval = 50;

	numberOfLakes = 10;	
	
	optimumPhenoMarine = 10;
	
	optimumPhenoFresh = -10;

	individualsInOcean = 2000;
	
	individualsInFresh = 2000;

	numberOfEffectRegions = 10;

	sizeOfEffectRegions = 100;

	sizeOfNonEffectBuffer = 4950;
	
	mutationRate = 1e-7;
		
	rrEffect = 1e-6;
	rrNeut = 1e-3;

	introduce = 40000;

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

	recombVector = NULL;
	lociVector = NULL;

	pos = 0;
	q = NULL;
	for (i in 1:C)
	{
		initializeGenomicElement(g1, pos, pos + W-1);    
		pos = pos + W;
		recombVector = c(recombVector,rrNeut);
		lociVector = c(lociVector,pos-1);

		initializeGenomicElement(g2, pos, pos + P-1);
		q = c(q,c(pos,pos+P-1));
		pos = pos + P;	
		recombVector = c(recombVector,rrEffect);
		lociVector = c(lociVector,pos-1);

		initializeGenomicElement(g1, pos, pos + W-1);
		pos = pos + W;
		recombVector = c(recombVector,rrNeut);
		lociVector = c(lociVector,pos-1);
	}

	catn("end of chromosome position: "+pos);
	totalNumOfLoci = pos;

	defineConstant("TL", totalNumOfLoci);
	defineConstant("Q", q);     // remember our QTL positions
		
	catn(Q);	
	catn(TL);

	initializeRecombinationRate(recombVector,lociVector);	

	//Interation Function for Freshwater mateChoice() callback...

	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.02);
	i1.setInteractionFunction("n", 0.5, 0.02);

	//Interaction Function for Oceanic mateChoice() callback...

	initializeInteractionType(2, "x", reciprocal=T, maxDistance=1.0);
	i2.setInteractionFunction("n", 0.5, 0.5);
	
	writeFile(PATH+"AveragePheno.txt",numberOfSamplingPoints+" "+INT+" "+((INP/INT)-1),append = F);
	//writeFile(PATH+"MeanVarience.txt",numberOfSamplingPoints+" "+INT,append = F);
	//writeFile(PATH+"AveragePheno.txt",numberOfSamplingPoints+" "+INT+" "+NLAKES,append = F);
	//writeFile(PATH+"myData.ped"," ",append = F);
	writeFile(PATH+"AverageFst.txt",numberOfSamplingPoints+" "+INT,append = F);
	writeFile(PATH+"AverageFst2.txt",numberOfSamplingPoints+" "+INT+" "+INP,append = F);
	writeFile(PATH+"OceanFreshwaterFst.txt"," ",append = F);
	writeFile(PATH+"OceanFreshwater2Fst.txt"," ",append = F);
	writeFile(PATH+"FreshwaterFreshwater2Fst.txt"," ",append = F);
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


40000 late(){

	sim.addSubpopSplit("p3",IIF,p1);
	p3.setSpatialBounds(c(0.0, LOP));
	
	lakeIndividuals2 = p3.individuals;
	
	p1.setMigrationRates(p3,MGR);
	p3.setMigrationRates(p1,MGR);

	for (i in 0:(NLAKES-1)){
		for (j in 0:(IPF-1)){
			lakeIndividuals2[IPF*i + j].setSpatialPosition((i+1)-(0.5));						
		}
	}
}
	

1: late() {  
	i1.evaluate();
	i2.evaluate();		
	
	substitutionOffset = 2 * sum(sim.substitutions.selectionCoeff);

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
	
		individual.tagF = dominantMarine+recessiveFreshwater+dominantFreshwater+recessiveMarine+additiveMarine+additiveFreshwater+substitutionOffset;		
	}


	if((sim.generation % INT == 0) | (sim.generation == END))
	{
		//TODO Number of mutations per individual. 
	
		//---------------------Sample Average Phenotype----------------------------------------
		cat("-------------------------------\n");
		cat("Output1 for generation " + sim.generation + ":\n\n");
		Oceanic_mean_phenotype =  mean(p1.individuals.tagF);
		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
		catn("Oceanic mean pheno: "+Oceanic_mean_phenotype);
		catn("Freshwater mean pheno: "+Freshwater_mean_phenotype);	


		if (sim.generation >= INP){
			newFreshwater_mean_phenotype = mean(p3.individuals.tagF);
			catn("new Freshwater mean pheno: "+newFreshwater_mean_phenotype);		
			writeFile(PATH+"AveragePheno.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype+" "+newFreshwater_mean_phenotype,append = T);
		}else{
			writeFile(PATH+"AveragePheno.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype,append = T);
		}


	}

	//Set tag equal to a value to be returned by fitness callback 
	
	for (individual in p2.individuals){
		individual.tagF = dnorm(-10.0 - individual.tagF,0.0,15.0);
	}
	for (individual in p1.individuals){
		individual.tagF = dnorm(10 - individual.tagF,0.0,15.0);
	}
	if (sim.generation >= INP){
		for (individual in p3.individuals){
			individual.tagF = dnorm(-10 - individual.tagF,0.0,15.0);
		}
	}

	//Normalize the fitness in the lakes so an equivelent mean number of parents get chosen from each.
 
	for (l in 1:NLAKES){
		locationOfLake = (l - 0.5);
		lake_individuals = p2.individuals[abs(p2.individuals.spatialPosition - locationOfLake) < 1e-8];
		sumOfPhenos = sum(lake_individuals.tagF);
		lake_individuals.tagF =(lake_individuals.tagF*IPF)/(sumOfPhenos);
	}
	
	if (sim.generation >= INP){
		for (l in 1:NLAKES){
			locationOfLake = (l - 0.5);
			lake_individuals2 = p3.individuals[abs(p3.individuals.spatialPosition - locationOfLake) < 1e-8];
			sumOfPhenos2 = sum(lake_individuals2.tagF);
			lake_individuals2.tagF =(lake_individuals2.tagF*IPF)/(sumOfPhenos2);
		}
	}

	if(sim.generation == END)
	{
		
		//FST PER MUTATION ----------------------------- OceanFreshwater ----------------------

		muts = sim.mutations;
		mutpos = muts.position;
		p1_p = sim.mutationFrequencies(p1,muts); //get frequencies of all mutations in terms of p1
		p2_p = sim.mutationFrequencies(p2,muts); //get frequencies of all mutations in terms of p2
		mean_p = (p1_p + p2_p) / 2.0;       //average the frequencies. 

		pq = (p1_p * (1 - p1_p) + p2_p * (1 - p2_p))/2.0;
		mean_pq = mean_p * (1 - mean_p);
		fst = 1 - (pq / mean_pq);
		mutpos = mutpos[!isNAN(fst)];
		fst = fst[!isNAN(fst)];

		writeFile(PATH+"OceanFreshwaterFst.txt",paste(fst),append = T);
		writeFile(PATH+"OceanFreshwaterFst.txt",paste(mutpos),append = T);

		//FST PER MUTATION ----------------------------- OceanFreshwater2 ----------------------

		muts = sim.mutations;
		mutpos = muts.position;
		p1_p = sim.mutationFrequencies(p1,muts); //get frequencies of all mutations in terms of p1
		p3_p = sim.mutationFrequencies(p3,muts); //get frequencies of all mutations in terms of p2
		mean_p = (p1_p + p3_p) / 2.0;       //average the frequencies. 

		pq = (p1_p * (1 - p1_p) + p3_p * (1 - p3_p))/2.0;
		mean_pq = mean_p * (1 - mean_p);
		fst = 1 - (pq / mean_pq);
		mutpos = mutpos[!isNAN(fst)];
		fst = fst[!isNAN(fst)];

		writeFile(PATH+"OceanFreshwater2Fst.txt",paste(fst),append = T);
		writeFile(PATH+"OceanFreshwater2Fst.txt",paste(mutpos),append = T);

		//FST PER MUTATION ----------------------------- FreshwaterFreshwater2 ----------------------

		muts = sim.mutations;
		mutpos = muts.position;
		p2_p = sim.mutationFrequencies(p2,muts); //get frequencies of all mutations in terms of p1
		p3_p = sim.mutationFrequencies(p3,muts); //get frequencies of all mutations in terms of p2
		mean_p = (p2_p + p3_p) / 2.0;       //average the frequencies. 

		pq = (p2_p * (1 - p2_p) + p3_p * (1 - p3_p))/2.0;
		mean_pq = mean_p * (1 - mean_p);
		fst = 1 - (pq / mean_pq);
		mutpos = mutpos[!isNAN(fst)];
		fst = fst[!isNAN(fst)];

		writeFile(PATH+"FreshwaterFreshwater2Fst.txt",paste(fst),append = T);
		writeFile(PATH+"FreshwaterFreshwater2Fst.txt",paste(mutpos),append = T);
		
		//-------------------Simulation Over---------------------

		sim.simulationFinished();
	}
}

function (f$)calcFST(o<Subpopulation>$ subpop1, o<Subpopulation>$ subpop2, o<Mutation> mutations)
{
	p1_p = sim.mutationFrequencies(subpop1,mutations);
	p2_p = sim.mutationFrequencies(subpop2,mutations);

	mean_p = (p1_p + p2_p) / 2.0;
	pq = (p1_p * (1 - p1_p) + p2_p * (1 - p2_p))/2.0;
	mean_pq = mean_p * (1 - mean_p);
	fst = 1 - (pq / mean_pq);		
	fst = fst[!isNAN(fst)];
	
	return mean(fst);
}

function (f$)calcFSTbetween(o<Subpopulation>$ subpop, o<Mutation> mutations)
{
	freqs = rep(0,size(mutations));
	mean_p = freqs;
	mean_het = freqs;

	for (i in 1:NLAKES)
	{	
		locationOfLake = (i - 0.5);
			
		lake_individuals = subpop.individuals[abs(subpop.individuals.spatialPosition - locationOfLake) < 1e-8];
		
		freqs = rep(0, size(mutations));
		
		for (geno in lake_individuals.genomes)
		{
			freqs = freqs + asInteger(geno.containsMutations(mutations));		
		}
		freqs = freqs / (2*size(lake_individuals));
		mean_p = mean_p + freqs;
		mean_het = mean_het + (freqs * (1-freqs));
	} 
	
	mean_p = mean_p/NLAKES;        
	mean_het = mean_het/NLAKES;         
	
	fst = (1 - mean_het / (mean_p * (1 - mean_p)));
	fst = mean(fst[!isNAN(fst)]);
	
	if(isNULL(fst)){
		fst = 0.0;
	}	
	
	return mean(fst);
}

function (f$)calcFSTbetweenMarine(o<Mutation> mutations)
{
	freqs = rep(0,size(mutations));
	mean_p = freqs;
	mean_het = freqs;

	for (i in 1:NLAKES)
	{	
		leftBound = (i - 1);
		rightBound = i;

		region_individuals = p1.individuals[p1.individuals.x >= leftBound & p1.individuals.x < rightBound];
					
		freqs = rep(0, size(mutations));
		
		for (geno in region_individuals.genomes)
		{
			freqs = freqs + asInteger(geno.containsMutations(mutations));		
		}
		freqs = freqs / (2*size(region_individuals));
		mean_p = mean_p + freqs;
		mean_het = mean_het + (freqs * (1-freqs));
	} 
	
	mean_p = mean_p/NLAKES;        
	mean_het = mean_het/NLAKES;         

	fst = (1 - mean_het / (mean_p * (1 - mean_p)));
	fst = mean(fst[!isNAN(fst)]);
	
	if(isNULL(fst)){
		fst = 0.0;
	}	

	return fst;
}

1: late(){		//FST SAMPLING P1 & P2


	if((sim.generation % INT == 0) | (sim.generation == END))
	{
		//GET MUTATIONS THAT ONLY EXIST IN WANTED POPULATIONS FOR SAMPLING.
		
		AverageFsts = NULL;

		m1s = sim.mutationsOfType(m1);			
		m2s = sim.mutationsOfType(m2);			
		m4s = sim.mutationsOfType(m4);			
		m5s = sim.mutationsOfType(m5);			
		m6s = sim.mutationsOfType(m6);			
		m7s = sim.mutationsOfType(m7);			
		
		allMutations = sim.mutations;
		effectMutations = c(m1s,m2s,m4s,m5s,m6s,m7s);
		neutralMutations = sim.mutationsOfType(m3);
		
		OceanFreshFst = calcFST(p1,p2,allMutations);
		AverageFsts = c(AverageFsts,OceanFreshFst);

		OceanFreshFstNeutral = calcFST(p1,p2,neutralMutations);
		AverageFsts = c(AverageFsts,OceanFreshFstNeutral);
		
		OceanFreshFstEffect = calcFST(p1,p2,effectMutations);
		AverageFsts = c(AverageFsts,OceanFreshFstEffect);
			 
		betweenLakesFst = calcFSTbetween(p2,allMutations);
		AverageFsts = c(AverageFsts,betweenLakesFst);

		betweenLakesFstNeutral = calcFSTbetween(p2,neutralMutations);	
		AverageFsts = c(AverageFsts,betweenLakesFstNeutral);

		betweenLakesFstEffect = calcFSTbetween(p2,effectMutations);	
		AverageFsts = c(AverageFsts,betweenLakesFstEffect);

		betweenMarineFst = calcFSTbetweenMarine(allMutations);	
		AverageFsts = c(AverageFsts,betweenMarineFst);

		betweenMarineFstNeutral = calcFSTbetweenMarine(neutralMutations);	
		AverageFsts = c(AverageFsts,betweenMarineFstNeutral);
	
		betweenMarineFstEffect = calcFSTbetweenMarine(effectMutations);	
		AverageFsts = c(AverageFsts,betweenMarineFstEffect);

		catn(AverageFsts);			
		writeFile(PATH+"AverageFst.txt",paste(AverageFsts),append = T);
	}			
}

40000: late(){

	if((sim.generation % INT == 0) | (sim.generation == END))
	{

		AverageFsts = NULL;

		m1s = sim.mutationsOfType(m1);			
		m2s = sim.mutationsOfType(m2);			
		m4s = sim.mutationsOfType(m4);			
		m5s = sim.mutationsOfType(m5);			
		m6s = sim.mutationsOfType(m6);			
		m7s = sim.mutationsOfType(m7);			

		allMutations = sim.mutations;
		effectMutations = c(m1s,m2s,m4s,m5s,m6s,m7s);
		neutralMutations = sim.mutationsOfType(m3);
	
		MarineIntroFst = calcFST(p1,p3,allMutations);	
		AverageFsts = c(AverageFsts,MarineIntroFst);
		
		MarineIntroFstNeutral = calcFST(p1,p3,neutralMutations);
		AverageFsts = c(AverageFsts,MarineIntroFstNeutral);
	
		MarineIntroFstEffect = calcFST(p1,p3,effectMutations);
		AverageFsts = c(AverageFsts,MarineIntroFstEffect);

		betweenIntroFst = calcFSTbetween(p3,allMutations);
		AverageFsts = c(AverageFsts,betweenIntroFst); 	

		betweenIntroFstNeutral = calcFSTbetween(p3,neutralMutations);
		AverageFsts = c(AverageFsts,betweenIntroFstNeutral); 	

		betweenIntroFstEffect = calcFSTbetween(p3,effectMutations);
		AverageFsts = c(AverageFsts,betweenIntroFstEffect); 	

		freshIntroFst = calcFST(p2,p3,allMutations);
		AverageFsts = c(AverageFsts,freshIntroFst);
		
		freshIntroFstNeutral = calcFST(p2,p3,neutralMutations);
		AverageFsts = c(AverageFsts,freshIntroFstNeutral);

		freshIntroFstEffect = calcFST(p2,p3,effectMutations);
		AverageFsts = c(AverageFsts,freshIntroFstEffect);

		writeFile(PATH+"AverageFst2.txt",paste(AverageFsts),append = T);
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
40000: fitness(NULL, p3){
	return individual.tagF;
}
1: mateChoice(p2) {
	return  i1.strength(individual) * weights; 
}
1: mateChoice(p1){
	return i2.strength(individual) * weights;
}
40000: mateChoice(p3){
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
40000: modifyChild(p3){

	if(sourceSubpop != subpop){  //migrant from Lake -> Oceanic 
		pos = parent1.spatialPosition + rnorm(1,0,0.5);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}else{
		child.setSpatialPosition(parent1.spatialPosition);
	}
	return T;
}



160001 late() {sim.simulationFinished();}
