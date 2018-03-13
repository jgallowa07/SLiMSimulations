//This model has mutations for both freshwater and marine that are recessive (complete) and dominant (complete)

initialize() {

	pathToOutputFile = "/projects/haldane/shared/SLiMSimulations/Output1/MyRecipe8_2_4/";
	//pathToOutputFile = "~/Documents/ResearchStuff/CODEBASS/Assigned/Output1/MyRecipe8_2_4Local/";
	//pathToOutputFile = "~/Desktop/";

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

	migrationRate = 0.0005;

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
		
//	catn(Q);	
//	catn(TL);

	initializeRecombinationRate(recombVector,lociVector);	

	//Interation Function for Freshwater mateChoice() callback...

	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.02);
	i1.setInteractionFunction("n", 0.5, 0.02);

	//Interaction Function for Oceanic mateChoice() callback...

	initializeInteractionType(2, "x", reciprocal=T, maxDistance=1.0);
	i2.setInteractionFunction("n", 0.5, 0.5);
	
	writeFile(PATH+"AveragePheno.txt",numberOfSamplingPoints+" "+INT+" "+((INP/INT)-1),append = F);
	writeFile(PATH+"EffectMutationsThrough.txt",numberOfSamplingPoints+" "+INT+" "+((INP/INT)-1),append=F);
	writeFile(PATH+"sharedAlleles.txt",numberOfSamplingPoints+" "+INT+" "+((INP/INT)-1),append=F);
	writeFile(PATH+"numfreshAlleles.txt",numberOfSamplingPoints+" "+INT+" "+((INP/INT)-1),append=F);
	writeFile(PATH+"AvgFreshAllelesPerMarineInd.txt",numberOfSamplingPoints+" "+INT+" "+((INP/INT)-1),append=F);
	writeFile(PATH+"OceanFreshwaterFst.txt"," ",append = F);
	writeFile(PATH+"OceanFreshwater2Fst.txt"," ",append = F);
	writeFile(PATH+"FreshwaterFreshwater2Fst.txt"," ",append = F);
	writeFile(PATH+"sgv.txt",(((END-INP)/5)+1)+" "+5+" "+INP+" ",append = F);
	writeFile(PATH+"avgFWAA.txt",(((END-INP)/5)+1)+" "+5+" "+INP+" ",append = F);
	
}

//--INTRODUCE MARINE AND FRESHWATER POPULATION--

1 late(){

	sim.setValue("Adapted",F);
	
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

	if (sim.generation >= INP & sim.generation % 5 == 0){

		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
		newFreshwater_mean_phenotype = mean(p3.individuals.tagF);
		
		sim.setValue("FAP",Freshwater_mean_phenotype);
		sim.setValue("IAP",newFreshwater_mean_phenotype);

	}

	//--SAMPLE PHENOTYPE--
	
	if((sim.generation % INT == 0) | (sim.generation == END))
	{
		//TODO Number of mutations per individual. 
	
//		cat("-------------------------------\n");
//		cat("Output1 for generation " + sim.generation + ":\n\n");
		Oceanic_mean_phenotype =  mean(p1.individuals.tagF);
		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
//		catn("Oceanic mean pheno: "+Oceanic_mean_phenotype);
//		catn("Freshwater mean pheno: "+Freshwater_mean_phenotype);	


		if (sim.generation >= INP){
			newFreshwater_mean_phenotype = mean(p3.individuals.tagF);
//			catn("new Freshwater mean pheno: "+newFreshwater_mean_phenotype);		
			writeFile(PATH+"AveragePheno.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype+" "+newFreshwater_mean_phenotype,append = T);
		}else{
			writeFile(PATH+"AveragePheno.txt",Oceanic_mean_phenotype+" "+Freshwater_mean_phenotype,append = T);
		}


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
		endsim();
	}
}

//--AT THE END OF THE SIM, SAMPLE FST ACROSS THE GENOME FOR ALL THREE COMBINATIONS OF SUBPOPS--

function (void)endsim(void){


		muts = sim.mutations;
		mutpos = muts.position;
		ids = muts.id;
		
		//FST PER MUTATION ----------------------------- OceanFreshwater ----------------------


		fst1 = calcFST(p1,p2,muts);

		ids1 = ids[!isNAN(fst1)];
		mutPosMatch1 = mutpos[!isNAN(fst1)];
		fst1 = fst1[!isNAN(fst1)];

		top_mutations1 = ids1[fst1 > 0.5];
		catn("{ Effect Mutations | fst(Ocean,Freshwater) > 0.5 } = " + paste(top_mutations1));

		writeFile(PATH+"OceanFreshwaterFst.txt",paste(fst1),append = T);
		writeFile(PATH+"OceanFreshwaterFst.txt",paste(mutPosMatch1),append = T);

		//FST PER MUTATION ----------------------------- OceanFreshwater2 ----------------------

		fst2 = calcFST(p1,p3,muts);

		ids2 = ids[!isNAN(fst2)];
		mutPosMatch2 = mutpos[!isNAN(fst2)];
		fst2 = fst2[!isNAN(fst2)];

		top_mutations2 = ids2[fst2 > 0.5];
		catn("{ Effect Mutations | fst(Ocean,Introduced) > 0.5 } = " + paste(top_mutations2));

		writeFile(PATH+"OceanFreshwater2Fst.txt",paste(fst2),append = T);
		writeFile(PATH+"OceanFreshwater2Fst.txt",paste(mutPosMatch2),append = T);

		//FST PER MUTATION ----------------------------- FreshwaterFreshwater2 ----------------------
		
		fst3 = calcFST(p2,p3,muts);

		mutPosMatch3 = mutpos[!isNAN(fst3)];
		fst3 = fst3[!isNAN(fst3)];

		writeFile(PATH+"FreshwaterFreshwater2Fst.txt",paste(fst3),append = T);
		writeFile(PATH+"FreshwaterFreshwater2Fst.txt",paste(mutPosMatch3),append = T);
	
		//------------------Effect Mutation Sampling-------------


		sim.simulationFinished();
	}


//return the frequecies for a set of mutations among a given set of individuals
//currently slim doesn't support this but I think it should ..
function (f)mutFreq(o<Mutation> mutations, o<Individual> individuals)
{
	freqs = rep(0,size(mutations));
	for (geno in individuals.genomes)
	{
		freqs = freqs + asInteger(geno.containsMutations(mutations));		
	}
	freqs = freqs / (2 * size(individuals));
	
	return freqs;
}

//--HELPER TO CALCULATE FST BETWEEN TWO SUBPOPULATIONS--

function (f)calcFST(o<Subpopulation>$ subpop1, o<Subpopulation>$ subpop2, o<Mutation> mutations)
{
	p1_p = sim.mutationFrequencies(subpop1,mutations);
	p2_p = sim.mutationFrequencies(subpop2,mutations);

	mean_p = (p1_p + p2_p) / 2.0;
	pq = (p1_p * (1 - p1_p) + p2_p * (1 - p2_p))/2.0;
	mean_pq = mean_p * (1 - mean_p);
	fst = 1 - (pq / mean_pq);		
	
	return fst;
}

//--HELPER TO GET ALL THE HIGH FREQ ALLELES FOR EACH LAKE, THAT ARE NOT HIGH FREQ IN MARINE

function (o<Mutation>)getHighFreqFreshAlleles(o<Subpopulation>$ subpop, o<Mutation> mutations)
{
	//filter out any alleles that are high frequency in the marine as well

	frequencyInMarine = sim.mutationFrequencies(p1,mutations);
	mutations = mutations[frequencyInMarine < 0.5];

	ret = NULL;	
	for (i in 1:NLAKES){

		locationOfLake = (i - 0.5);
		lake_individuals = subpop.individuals[abs(subpop.individuals.spatialPosition - locationOfLake) < 1e-8];
		freq = mutFreq(mutations,lake_individuals);
		Ki = mutations[freq > 0.5];
		ret = c(ret,Ki);		

	}	
	
	return ret;
}

//--RECORD THE EFFECT MUTATION FREQUENCIES IN ALL POPS--

function (void)recordFrequencies(o<Mutation> mutations){

	if(sim.generation == INP){
		writeFile(PATH+"EffectMutationsThrough.txt","#",append = T);
	}

	Ids = mutations.id;
	writeFile(PATH+"EffectMutationsThrough.txt",paste(Ids),append = T);	
	freq1 = sim.mutationFrequencies(p1,mutations);
	freq2 = sim.mutationFrequencies(p2,mutations);
	writeFile(PATH+"EffectMutationsThrough.txt",paste(freq1),append = T);
	writeFile(PATH+"EffectMutationsThrough.txt",paste(freq2),append = T);

	if(sim.generation >= INP){
		freq3 = sim.mutationFrequencies(p3,mutations);
		writeFile(PATH+"EffectMutationsThrough.txt",paste(freq3),append = T);
	}
	
	return;

}

//--RECORD 1. AVG NUMBER OF LAKES EACH FWAA ALLELE APPEARS IN 2. AVG NUMBER OF FWAA / MARINE IND 3. NUM FWAA 

function (void)recordFreshAdaptedAlleles(o<Mutation> mutations){

	highFreqAlleles = getHighFreqFreshAlleles(p2,mutations);
	freshwaterAdaptedAlleles = unique(highFreqAlleles);
	avgSharedAlleles = size(highFreqAlleles) / size(freshwaterAdaptedAlleles);
	avgSharedAlleles[isNAN(avgSharedAlleles)] = 0;


	freqInMarine = sim.mutationFrequencies(p1,freshwaterAdaptedAlleles);	
	avgFreshAllelesPerMarineInd = sum(freqInMarine);

	writeFile(PATH+"sharedAlleles.txt",asString(avgSharedAlleles),append = T);
	writeFile(PATH+"AvgFreshAllelesPerMarineInd.txt",asString(avgFreshAllelesPerMarineInd),append = T);
	writeFile(PATH+"numfreshAlleles.txt",asString(size(freshwaterAdaptedAlleles)),append=T);

	return;

}

//--RECORD THE SGV IN A POPULATION

function (void)recordStandingGeneticVariation(o<Subpopulation> subpop, o<Mutation> mutations){

	sgv_s = NULL;

	for (sub in subpop){		

		numEffMuts = size(mutations);
		total = 0;

		for (ef in mutations){
			s = ef.selectionCoeff;
			p = sim.mutationFrequencies(sub,ef);
			additiveGeneticVariance = (s*s)*(p)*(1-p);
			total = total + additiveGeneticVariance;
		}

		sgv = total/numEffMuts;
		sgv_s = c(sgv_s,sgv);
	}

	writeFile(PATH+"sgv.txt",paste(sgv_s),append = T);
	return;

}

//--SAMPLING--

1: late(){		

	if((sim.generation % INT == 0) | (sim.generation == END))
	{
		//GET MUTATIONS THAT ONLY EXIST IN WANTED POPULATIONS FOR SAMPLING.
		
		AverageFsts = NULL;
		highFreqAllelesInp2 = NULL;
		highFreqAllelesInp3 = NULL;

		m1s = sim.mutationsOfType(m1);			
		m2s = sim.mutationsOfType(m2);			
		m4s = sim.mutationsOfType(m4);			
		m5s = sim.mutationsOfType(m5);			
		m6s = sim.mutationsOfType(m6);			
		m7s = sim.mutationsOfType(m7);			
		
		effectMutations = c(m1s,m2s,m4s,m5s,m6s,m7s);

		recordFrequencies(effectMutations);
		
		recordFreshAdaptedAlleles(effectMutations);
		
	}			
}

//--THIS IS WHERE WE HAVE A MORE VERBODE SAMPLING INTERVAL, 5, AND TRY TO ANALYZE THE RECONSTRUCTION OF THE FRESHWATER GENOTYPE 

40000: late(){
	
	if (sim.generation % 5 == 0){
		
		m1s = sim.mutationsOfType(m1);			
		m2s = sim.mutationsOfType(m2);			
		m4s = sim.mutationsOfType(m4);			
		m5s = sim.mutationsOfType(m5);			
		m6s = sim.mutationsOfType(m6);			
		m7s = sim.mutationsOfType(m7);			
		
		effectMutations = c(m1s,m2s,m4s,m5s,m6s,m7s);
				
		origFreshMeanPheno = sim.getValue("FAP");
		introFreshMeanPheno = sim.getValue("IAP");

//		catn(sim.generation + " " + origFreshMeanPheno);
//		catn(sim.generation + " " + introFreshMeanPheno);

		highFreqAlleles = getHighFreqFreshAlleles(p2,effectMutations);
		FWAA = unique(highFreqAlleles);

		avgFWAAinMarine = sum(sim.mutationFrequencies(p1,FWAA));
		avgFWAAinOrigFresh = sum(sim.mutationFrequencies(p2,FWAA));
		avgFWAAinIntroFresh = sum(sim.mutationFrequencies(p3,FWAA));

		writeFile(PATH+"avgFWAA.txt",paste(c(avgFWAAinMarine,avgFWAAinOrigFresh,avgFWAAinIntroFresh)),append = T);		

		recordStandingGeneticVariation(c(p1,p2,p3),effectMutations);

		if (abs(origFreshMeanPheno - introFreshMeanPheno) <= 1 & !(sim.getValue("Adapted"))){
			
			timeToAdaptation = sim.generation - INP;			
			writeFile(PATH+"timeToAdaptation.txt",timeToAdaptation + " ",append = F);
			
			highFreqAllelesInIntro = getHighFreqFreshAlleles(p3,effectMutations);
			FWAAIN = unique(highFreqAllelesInIntro);

			p2_ids = FWAA.id;
			p3_ids = FWAAIN.id;

			writeFile(PATH+"FWAA_IDs_ORIG_INTRO",paste(p2_ids),append = F);
			writeFile(PATH+"FWAA_IDs_ORIG_INTRO",paste(p3_ids),append = T);

			sim.setValue("Adapted",T);
		}
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
40000: fitness(NULL, p3){
	return individual.tagF;
}

//--CHOOSE MATES BASED UPON FITNESS & SPATIAL LOCATION--

1: mateChoice(p2) {
	return  i1.strength(individual) * weights; 
}
1: mateChoice(p1){
	return i2.strength(individual) * weights;
}
40000: mateChoice(p3){
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
