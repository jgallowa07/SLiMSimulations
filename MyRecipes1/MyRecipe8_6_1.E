//This model has mutations for both freshwater and marine that are recessive (complete) and dominant (complete)

initialize() {

	//Source File With Helper Functions
	//source("Utils.slim");
	source("/projects/haldane/shared/SLiMSimulations/MyRecipes1/Utils.slim");

	//The Path Where You want sampling files written to
	pathToOutputFile = "/projects/haldane/shared/SLiMSimulations/Output1/MyRecipe8_6_1/";
	//pathToOutputFile = "~/Desktop/";

	//KNOBS//

	endOfSimulation = 80000;

	introduce = 40000;

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
	
	recombinationRateInEffectRegions = 1e-5;

	recombinationRateInNeutralRegions = 1e-5;

	migrationRate = 5e-04;

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
		recombVector = c(recombVector,recombinationRateInNeutralRegions);
		lociVector = c(lociVector,pos-1);

		initializeGenomicElement(g2, pos, pos + P-1);
		q = c(q,c(pos,pos+P-1));
		pos = pos + P;	
		recombVector = c(recombVector,recombinationRateInEffectRegions);
		lociVector = c(lociVector,pos-1);

		initializeGenomicElement(g1, pos, pos + W-1);
		pos = pos + W;
		recombVector = c(recombVector,recombinationRateInNeutralRegions);
		lociVector = c(lociVector,pos-1);
	}

	catn("end of chromosome position: "+pos);
	totalNumOfLoci = pos;

	defineConstant("TL", totalNumOfLoci);
	defineConstant("Q", q);     // remember our QTL positions

	initializeRecombinationRate(recombVector,lociVector);	

	//Interation Function for Freshwater mateChoice() callback...

	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.02);
	i1.setInteractionFunction("n", 0.5, 0.02);

	//Interaction Function for Oceanic mateChoice() callback...

	initializeInteractionType(2, "x", reciprocal=T, maxDistance=1.0);
	i2.setInteractionFunction("n", 0.5, 0.5);


	//file header for ALL sampling files = numberOfSamplingPoints + sampingInterval + introduction + endofsim
	header = numberOfSamplingPoints+" "+INT+" "+INP+" "+END;
	
	writeFile(PATH+"AveragePhenotypeThroughout.txt",header,append = F);
//	writeFile(PATH+"EffectMutationFrequenciesThroughout.txt",header,append=F);	
	writeFile(PATH+"AvgNumLakesFreshAlleleAppearsIn.txt",header,append=F);
	writeFile(PATH+"numfreshAlleles.txt",header,append=F);
	writeFile(PATH+"AvgFreshAllelesPerMarineInd.txt",header,append=F);
	writeFile(PATH+"OceanFreshwaterFst.txt",header,append = F);
	writeFile(PATH+"OceanFreshwater2Fst.txt",header,append = F);
	writeFile(PATH+"FreshwaterFreshwater2Fst.txt",header,append = F);
	writeFile(PATH+"StandingGeneticVariationThroughout.txt",header,append = F);
	writeFile(PATH+"avgFWAA.txt",header,append = F);	
	writeFile(PATH+"avgFWAA_divTotal.txt",header,append = F);	
	writeFile(PATH+"Adaptation.txt",header,append = F);	
	writeFile(PATH+"MeanFstThroughoutNeutral.txt",header,append = F);
	writeFile(PATH+"MeanFstThroughoutEffect.txt",header,append = F);
	writeFile(PATH+"HaplotypeVectors",header,append = F);
	writeFile(PATH+"FWAAPerMarineIndCounts",header,append = F);

	//TODO CHECK ALL FILE PATHS
	
}

//right off the bat, set adapted to false because the new populations have not adapted
//also schedule the introduction blocks to the defined value in init
//--CUT MIGRATION RIGHT BEFORE INTRODUCTION TO ENSURE NO IMMEDIATE MIGRANTS
39999 late(){
	p1.setMigrationRates(p2,0.0);
	p2.setMigrationRates(p1,0.0);
}
//--BRING MIGRATION BACK AFTER
40000 late(){
	p1.setMigrationRates(p2,MGR);
	p2.setMigrationRates(p1,MGR);
}

1 early(){
	sim.setValue("Adapted",F);

	sim.rescheduleScriptBlock(s1,INP,INP);

	for (block in c(s2,s3,s4,s5)){
		sim.rescheduleScriptBlock(block,INP);
	}	
}

//--INTRODUCE MARINE AND FRESHWATER POPULATION--

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

//--INTRODUCE NEW FRESH POPULATION--

s1 40000 late(){

	//before introduction, lets sample the ditribution of FWAA in Marine individuals.
	recordNumberOfFWAAPerMarineInd();

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

	//--EVALUATE SPATIAL INTERALCTION FUNCTIONS--
	
	i1.evaluate();
	i2.evaluate();		
	
	//--CALCULATE PHENOTYPE--
		
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

	//every generation, get reference to average phenotype.

	marine_mean_phenotype = mean(p1.individuals.tagF);
	freshwater_mean_phenotype = mean(p2.individuals.tagF);

	sim.setValue("FAP",freshwater_mean_phenotype);
	sim.setValue("MAP",marine_mean_phenotype);

	if (sim.generation >= INP){
		introduced_mean_phenotype = mean(p3.individuals.tagF);
		sim.setValue("IAP",introduced_mean_phenotype);
	}
	
	//--SAMPLE PHENOTYPE--
	
	if((sim.generation % INT == 0) | (sim.generation == END))
	{
		//TODO Number of mutations per individual. 
	
		cat("-------------------------------\n");
		cat("Output1 for generation " + sim.generation + ":\n\n");

	
		samplePhenotype();
			
	}

	//--CALCULATE FITNESS BASED UPON PHENOTYPE--
	
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

	//--NORMALIZE FITNESS IN LAKES SO PROPOTIONAL AMOUNT GET CHOSEN IN THE NEXT GEN--
 
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

	//--SAMPLE FST ACROSS THE GENOME AT THE END BETWEEN ALL PAIRS OF SUBPOPS--

	if(sim.generation == END)
	{
		FstAcrossChromosome();
		sim.simulationFinished();
	}
}

//--AT THE END OF THE SIM, SAMPLE FST ACROSS THE GENOME FOR ALL THREE COMBINATIONS OF SUBPOPS--

1: late(){		

	if((sim.generation % INT == 0) | (sim.generation == END))
	{
		//GET MUTATIONS THAT ONLY EXIST IN WANTED POPULATIONS FOR SAMPLING.

		m1s = sim.mutationsOfType(m1);			
		m2s = sim.mutationsOfType(m2);			
		m4s = sim.mutationsOfType(m4);			
		m5s = sim.mutationsOfType(m5);			
		m6s = sim.mutationsOfType(m6);			
		m7s = sim.mutationsOfType(m7);			
		
		effectMutations = c(m1s,m2s,m4s,m5s,m6s,m7s);

		//recordFrequencies(effectMutations);
		
		recordFreshAdaptedAlleles(effectMutations);
		
	}			
}

//--ANALYZE THE RECONSTRUCTION OF THE FRESHWATER GENOTYPE 

s2 40000: late(){

	//TODO RECORD P2 HAPLO AT TIME OF INTRODUCTION ... 


	//time of adaptation
	if (abs(sim.getValue("FAP") - sim.getValue("IAP")) <= 0.5 & !(sim.getValue("Adapted"))){

		sim.setValue("Adapted",T);
		
		//record time it took to adapt
		timeToAdaptation = sim.generation - INP;			

		
		m1s = sim.mutationsOfType(m1);			
		m2s = sim.mutationsOfType(m2);			
		m4s = sim.mutationsOfType(m4);			
		m5s = sim.mutationsOfType(m5);			
		m6s = sim.mutationsOfType(m6);			
		m7s = sim.mutationsOfType(m7);			

		allMutations = sim.mutations;
		effectMutations = c(m1s,m2s,m4s,m5s,m6s,m7s);
		neutralMutations = sim.mutationsOfType(m3);		
		
		//get high frequency alleles in original lakes
		highFreqAllelesInOrig = getHighFreqFreshAlleles(p2,effectMutations);
		FWAA_P2 = unique(highFreqAllelesInOrig);

		//get high frequency alleles in introduced lakes
		highFreqAllelesInIntro = getHighFreqFreshAlleles(p3,effectMutations);
		FWAA_P3 = unique(highFreqAllelesInIntro);
		
		numberOfSharedFWAA = sum(match(FWAA_P2.id,FWAA_P3.id) >= 0);

		totalNumOfFWAA = size(unique(c(FWAA_P2,FWAA_P3)));

		freqFWAA_P2InP2 = sim.mutationFrequencies(p2,FWAA_P2);
		freqFWAA_P2InP3 = sim.mutationFrequencies(p3,FWAA_P2);
		
		corFWAA_P2_freqs = cor(freqFWAA_P2InP2,freqFWAA_P2InP3);	
		
		freqEffectInP2 = sim.mutationFrequencies(p2,effectMutations);	
		freqEffectInP3 = sim.mutationFrequencies(p3,effectMutations);

		corEffect_freqs = cor(freqEffectInP2,freqEffectInP3);
			
		writeFile(PATH+"Adaptation.txt",paste(FWAA_P2.id),append = T);
		writeFile(PATH+"Adaptation.txt",paste(FWAA_P3.id),append = T);
		writeFile(PATH+"Adaptation.txt",asString(timeToAdaptation),append = T);
		writeFile(PATH+"Adaptation.txt",asString(numberOfSharedFWAA),append = T);
		writeFile(PATH+"Adaptation.txt",asString(totalNumOfFWAA),append = T);
		writeFile(PATH+"Adaptation.txt",asString(corFWAA_P2_freqs),append = T);
		writeFile(PATH+"Adaptation.txt",asString(corEffect_freqs),append = T);
		
		sim.outputFull(PATH+"simStateAtAdaptation.txt");
		

		// TODO RECORD FREQUENCIES AND CORRELATIONS ....

	}	

	if((sim.generation % INT == 0) | (sim.generation == END)){
		
		m1s = sim.mutationsOfType(m1);			
		m2s = sim.mutationsOfType(m2);			
		m4s = sim.mutationsOfType(m4);			
		m5s = sim.mutationsOfType(m5);			
		m6s = sim.mutationsOfType(m6);			
		m7s = sim.mutationsOfType(m7);			

		effectMutations = c(m1s,m2s,m4s,m5s,m6s,m7s);

		highFreqAlleles = getHighFreqFreshAlleles(p2,effectMutations);
		FWAA_P2 = unique(highFreqAlleles);

		numFWAA_P2 = size(FWAA_P2);

		avgFWAAinMarine = sum(sim.mutationFrequencies(p1,FWAA_P2));
		avgFWAAinOrigFresh = sum(sim.mutationFrequencies(p2,FWAA_P2));
		avgFWAAinIntroFresh = sum(sim.mutationFrequencies(p3,FWAA_P2));

		writeFile(PATH+"avgFWAA.txt",paste(c(avgFWAAinMarine,avgFWAAinOrigFresh,avgFWAAinIntroFresh)),append = T);

		PercentFWAAinMarine = avgFWAAinMarine / numFWAA_P2;
		PercentFWAAinOrigFresh = avgFWAAinOrigFresh / numFWAA_P2;
		PercentFWAAinIntroFresh = avgFWAAinIntroFresh / numFWAA_P2;

		if (!numFWAA_P2 == 0){	
			writeFile(PATH+"avgFWAA_divTotal.txt",paste(c(PercentFWAAinMarine,PercentFWAAinOrigFresh,PercentFWAAinIntroFresh)),append = T);	
		}else{
			writeFile(PATH+"avgFWAA_divTotal.txt",paste(c(0,0,0)),append = T);
		}

		recordFsts();
		recordStandingGeneticVariation(c(p1,p2,p3),effectMutations);	
	}
	
}




//--TO OVERRIDE SLIM'S DEFAULT DISTRIBUTION OF FITNESS EFFECTS--

fitness(m1) {return 1.0;}
fitness(m2) {return 1.0;}
fitness(m4) {return 1.0;}
fitness(m5) {return 1.0;}
fitness(m6) {return 1.0;}
fitness(m7) {return 1.0;}

//--FITNESS == 
//--WE DEFINE FITNESS TO BE DIRECTLY CORRELATED TO THIER PHENOTYPE--

fitness(NULL, p1) {
	return individual.tagF;
}
fitness(NULL, p2) {
	return individual.tagF;
}
s3 40000: fitness(NULL, p3){
	return individual.tagF;
}

//--CHOOSE MATES BASED UPON FITNESS & SPATIAL LOCATION--

1: mateChoice(p2) {
	return  i1.strength(individual) * weights; 
}
1: mateChoice(p1){
	return i2.strength(individual) * weights;
}
s4 40000: mateChoice(p3){
	return i1.strength(individual) * weights;
}

//--IF (THE OFFSPRING IS A MIGRANT FROM THE OCEAN) 
       	//FIND THE CLOSEST LAKE TO MIGRATE TO
//--ELSE 
	//MIGRATION FROM PARENT PULLED FROM NORMAL DISTRIBUTION
 
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

//--IF(THE OFFSPRING IS MIGRANT FROM THE LAKE)
	//MIGRATION FROM PARENT'S LAKE PULLED FROM NORMAL DIASTRIBUTION
//--ELSE
	//OFFSPRING IS PLACED IN EXACTLY THE SAME LAKE (POSITION) AS PARENT

modifyChild(p2){

	if(sourceSubpop != subpop){  //migrant from Lake -> Oceanic 
		pos = parent1.spatialPosition + rnorm(1,0,0.5);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}else{
		child.setSpatialPosition(parent1.spatialPosition);
	}
	return T;
}

s5 40000: modifyChild(p3){

	if(sourceSubpop != subpop){  //migrant from Lake -> Oceanic 
		pos = parent1.spatialPosition + rnorm(1,0,0.5);
		child.setSpatialPosition(p1.pointReflected(pos));	
	}else{
		child.setSpatialPosition(parent1.spatialPosition);
	}
	return T;
}

160001 late() {sim.simulationFinished();}
