initialize() {


	endOfSimulation = 10000;
	samplingInterval = 100;

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

	
	initializeMutationType("m1", 0.5, "g", 0.5, 1.0);  //Oceanic Benefitial Mutations
	initializeMutationType("m2", 0.5, "g", -0.5, 1.0); //Freshwater Benefictial Mutations
	initializeMutationType("m3", 0.5, "f", 0.0);       //Nuetral Benefitial Mutations
	initializeMutationType("m4", 0.5, "g", 0.0, 1.0);  //Contribution Mutations 
	
	initializeGenomicElementType("g1", m3, 1.0);
	initializeGenomicElementType("g2", c(m1,m2), c(1.0,1.0));
	initializeGenomicElementType("g3", m4, 1.0);

//	initializeGenomicElement(g2,0,999);

	m1.mutationStackPolicy = "l";
	m2.mutationStackPolicy = "l";

	defineConstant("C", 5);    // number of QTLs
	defineConstant("W", 100);  // size of neutral buffer on each side
//	defineConstant("P", 20);   //Size of contribution (Padding) Loci 

	
	//initialize genomic arch.
	// |----------Nuetral-----------|-Contributing-|QTL|-Contributing-|----------Nuetral----------|	  Repeat for number of QTL's ...



	pos = 0;
	q = NULL;
	
	for (i in 1:C)
	{
		initializeGenomicElement(g1, pos, pos + W-1);    
		pos = pos + W;
	
//		initializeGenomicElement(g3, pos, pos + P-1);
//		pos = pos + P;	
	
		initializeGenomicElement(g2, pos, pos);
		q = c(q, pos);
		pos = pos + 1;

//		initializeGenomicElement(g3, pos, pos + P-1);
//		pos = pos + P;	
		
		initializeGenomicElement(g1, pos, pos + W-1);
		pos = pos + W;
	}
	
	defineConstant("Q", q);     // remember our QTL positions
//	rates = c(rep(c(1e-8, 0.5), C-1), 1e-8);
//	ends = (repEach(Q + W, 2) + rep(c(0,1), C))[0:(C*2 - 2)];
	initializeRecombinationRate(1e-3);	

	//Interation Function for mateChoice() callback...

	initializeInteractionType(1, "x", reciprocal=T, maxDistance=0.02);
	i1.setInteractionFunction("n", 0.5, 0.02);
	
	
	writeFile("./myOut197_Pheno.txt",numberOfSamplingPoints+" "+samplingInterval,append = F);
	
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

}
	

1 late() {

	//choose 100 locations along the chromosome, give each genome a mutation of type m1 or m2 at that location
 
	g = sim.subpopulations.genomes;
	n = size(g);
	//mutations  = asInteger(round(runif(10,0,999)));
	//mutationSites = sort(mutations);
	
	for (q in Q)
	{
		isPlus = asLogical(rbinom(n, 1, 0.5));
		g[isPlus].addNewMutation(m1, 1.0, q);
		g[!isPlus].addNewMutation(m2, -1.0, q);
	}
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
		EffectMuts = c(individual.genomes.mutationsOfType(m1),individual.genomes.mutationsOfType(m2));
		individual.tagF = substitutionOffset + sum(EffectMuts.selectionCoeff);		
	}
	

	//The entire if block below is for sampling data at a given interval, INT	
	
	if(((sim.generation-1) % INT == 0) | (sim.generation == END))
	{
		cat("-------------------------------\n");
		cat("Output for generation " + sim.generation + ":\n\n");
		Oceanic_mean_phenotype =  mean(p1.individuals.tagF);
		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
		catn("Oceanic mean pheno: "+Oceanic_mean_phenotype);
		catn("Freshwater mean pheno: "+Freshwater_mean_phenotype);			
	
		writeFile("./myOut197_Pheno.txt",Oceanic_mean_phenotype + " " + Freshwater_mean_phenotype,append = T);
		
	}
	if(sim.generation == END)
	{
		
		Oceanic_mean_phenotype =  mean(p1.individuals.tagF);
		Freshwater_mean_phenotype = mean(p2.individuals.tagF);
		catn("Oceanic mean pheno"+Oceanic_mean_phenotype);
		catn("Freshwater mean pheno"+Freshwater_mean_phenotype);			

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
		writeFile("./myOut197_Muts.txt",num+" ",append = F);
		for (i in 0:(num-1))
		{
			writeFile("./myOut197_Muts.txt",mutfreq[i]+" "+mutpos[i]+" "+muteff[i],append = T);
		}

		mutfreq1 = sim.mutationFrequencies(p1,p1muts);
		muteff1 = p1muts.selectionCoeff;
		mutpos1 = p1muts.position;
		writeFile("./myOut197_1Muts.txt",p1num+" ",append = F);
		for (i in 0:(p1num-1))
		{
			writeFile("./myOut197_1Muts.txt",mutfreq1[i]+" "+mutpos1[i]+" "+muteff1[i],append = T);
		}

		mutfreq2 = sim.mutationFrequencies(p2,p2muts);
		muteff2 = p2muts.selectionCoeff;
		mutpos2 = p2muts.position;
		writeFile("./myOut197_2Muts.txt",p2num+" ",append = F);
		for (i in 0:(p2num-1))
		{
			writeFile("./myOut197_2Muts.txt",mutfreq2[i]+" "+mutpos2[i]+" "+muteff2[i],append = T);
		}

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
