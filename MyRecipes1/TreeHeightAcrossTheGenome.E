
initialize() {


	//Source File With Helper Functions
	//source("Utils.slim");
	//source("/projects/haldane/shared/SLiMSimulations/MyRecipes1/Utils.slim");
	//source("/Users/jaredgalloway/Documents/ResearchStuff/CODEBASS/Assigned/MyRecipes1/Utils.slim")

	//The Path Where You want sampling files written to
	pathToOutputFile = "/projects/haldane/shared/SLiMSimulations/Output1/TreeHeightAcrossTheGenome/";
	//pathToOutputFile = "~/Documents/ResearchStuff/CODEBASS/Assigned/Output1/MyRecipe8_3_0Local/";
	pathToOutputFile = "~/Desktop/";

	setSeed(1527001974015);

	//KNOBS//

	endOfSimulation = 80;

	introduce = 100000;

	samplingInterval = 100000;

	numberOfLakes = 10;	
	
	optimumPhenoMarine = 10;
	
	optimumPhenoFresh = -10;

	individualsInOcean = 200;
	
	individualsInFresh = 200;

	numberOfEffectRegions = 10;

	sizeOfEffectRegions = 100;

	sizeOfNonEffectBuffer = 4950;
	
	mutationRateInEffectRegions = 1e-7;
	
	mutationRateInNeutralRegions = 0;
		
	recombinationRateInEffectRegions = 1e-3;

	recombinationRateInNeutralRegions = 1e-3;

	migrationRate = 5e-05;

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
	//defineConstant("MR", mutationRate);
	defineConstant("C", numberOfEffectRegions);    	// number of Effect Regions;
	defineConstant("W", sizeOfNonEffectBuffer);  	// size of neutral buffer on each side
	defineConstant("P", sizeOfEffectRegions);   	// Size of contribution (Effect) Loci 
	defineConstant("INP", introduce);
	defineConstant("MGR", migrationRate);

	initializeSLiMOptions(dimensionality="x");

	initializeTreeSeq();

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

	mutRateVector = NULL;
	recombVector = NULL;
	lociVector = NULL;

	pos = 0;
	q = NULL;
	for (i in 1:C)
	{
		initializeGenomicElement(g1, pos, pos + W-1);    
		pos = pos + W;
		recombVector = c(recombVector,recombinationRateInNeutralRegions);
		mutRateVector = c(mutRateVector,mutationRateInNeutralRegions);
		lociVector = c(lociVector,pos-1);

		initializeGenomicElement(g2, pos, pos + P-1);
		q = c(q,c(pos,pos+P-1));
		pos = pos + P;	
		recombVector = c(recombVector,recombinationRateInEffectRegions);
		mutRateVector = c(mutRateVector,mutationRateInEffectRegions);
		lociVector = c(lociVector,pos-1);

		initializeGenomicElement(g1, pos, pos + W-1);
		pos = pos + W;
		recombVector = c(recombVector,recombinationRateInNeutralRegions);
		mutRateVector = c(mutRateVector,mutationRateInNeutralRegions);
		lociVector = c(lociVector,pos-1);
	}

	catn("end of chromosome position: "+pos);
	totalNumOfLoci = pos;

	defineConstant("TL", totalNumOfLoci);
	defineConstant("Q", q);     // remember our QTL positions

	initializeRecombinationRate(recombVector,lociVector);	
	initializeMutationRate(mutRateVector,lociVector); 

	//Interation Function for Freshwater mateChoice() callback...

	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.02);
	i1.setInteractionFunction("n", 0.5, 0.02);

	//Interaction Function for Oceanic mateChoice() callback...

	initializeInteractionType(2, "x", reciprocal=T, maxDistance=1.0);
	i2.setInteractionFunction("n", 0.5, 0.5);


	//file header for ALL sampling files = numberOfSamplingPoints + sampingInterval + introduction + endofsim
	header = numberOfSamplingPoints+" "+INT+" "+INP+" "+END;
/*
	writeFile(PATH+"AveragePhenotypeThoughout.txt",header,append = F);
	writeFile(PATH+"EffectMutationFrequenciesThroughout.txt",numberOfSamplingPoints+" "+INT+" "+((INP/INT)-1),append=F);	
	writeFile(PATH+"AvgNumLakesFreshAlleleAppearsIn.txt",header,append=F);
	writeFile(PATH+"numfreshAlleles.txt",header,append=F);
	writeFile(PATH+"AvgFreshAllelesPerMarineInd.txt",header,append=F);
	writeFile(PATH+"OceanFreshwaterFst.txt",header,append = F);
	writeFile(PATH+"OceanFreshwater2Fst.txt",header,append = F);
	writeFile(PATH+"FreshwaterFreshwater2Fst.txt",header,append = F);
//	writeFile(PATH+"StandingGeneticVariationThroughout.txt",header,append = F);
	writeFile(PATH+"avgFWAA.txt",header,append = F);	
	writeFile(PATH+"avgFWAA_divTotal.txt",header,append = F);	
	writeFile(PATH+"Adaptation.txt",header,append = F);	
	writeFile(PATH+"MeanFstThroughoutNeutral.txt",header,append = F);
	writeFile(PATH+"MeanFstThroughoutEffect.txt",header,append = F);
*/
	//TODO CHECK ALL FILE PATHS
	
}

//right off the bat, set adapted to false because the new populations have not adapted
//also schedule the introduction blocks to the defined value in init

//--INTRODUCE MARINE AND FRESHWATER POPULATION--

1 late(){	
	
	sim.addSubpop("p1", IIO);  //Marine 
	
	//set length of pops equal to 1 unit per lake

	p1.setSpatialBounds(c(0.0, LOP));

	//set migration to 1 individual per generation from each pop to another 

	//unifomily distribute oceanic individuals.

	for (ind in p1.individuals){
		ind.setSpatialPosition(p1.pointUniform());
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

	//--SAMPLE PHENOTYPE--
/*	
	if((sim.generation % INT == 0) | (sim.generation == END))
	{
		//TODO Number of mutations per individual. 
	
//		cat("-------------------------------\n");
//		cat("Output1 for generation " + sim.generation + ":\n\n");
//		catn("Oceanic mean pheno: "+sim.getValue("MAP"));
//		catn("Freshwater mean pheno: "+sim.getValue("FAP"));	

		if (sim.generation >= INP){
//			catn("new Freshwater mean pheno: "+sim.getValue("IAP"));		
//			writeFile(PATH+"AveragePhenotypeThoughout.txt",sim.getValue("MAP")+" "+sim.getValue("FAP")+" "+sim.getValue("IAP"),append = T);
		}else{
//			writeFile(PATH+"AveragePhenotypeThoughout.txt",sim.getValue("MAP")+" "+sim.getValue("FAP"),append = T);
		}
	}
*/
	//--CALCULATE FITNESS BASED UPON PHENOTYPE--
	for (individual in p1.individuals){
		individual.tagF = dnorm(10 - individual.tagF,0.0,15.0);
	}
	//--SAMPLE FST ACROSS THE GENOME AT THE END BETWEEN ALL PAIRS OF SUBPOPS--

	if(sim.generation == END)
	{

		sim.treeSeqOutput(PATH,binary = F,simplify = T);
		sim.simulationFinished();
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
//--CHOOSE MATES BASED UPON FITNESS & SPATIAL LOCATION--

1: mateChoice(p1) {
	return  i1.strength(individual) * weights; 
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


160001 late() {sim.simulationFinished();}
