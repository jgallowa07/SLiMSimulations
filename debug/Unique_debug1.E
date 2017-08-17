initialize() {
	initializeMutationRate(1e-7);
	
	// m1 mutation type: neutral
	initializeMutationType("m1", 0.5, "f", 0.0);
	
	// g1 genomic element type: uses m1 for all mutations
	initializeGenomicElementType("g1", m1, 1.0);
	
	// uniform chromosome of length 100 kb with uniform recombination
	initializeGenomicElement(g1, 0, 99999);
	initializeRecombinationRate(1e-8);
}

// create a population of 500 individuals
1 {
	sim.addSubpop("p1", 1500);
}

1 late() {

	//choose 100 locations along the chromosome, give each genome a mutation of type m1 or m2 at that location
 
	G = sim.subpopulations.genomes;
	n = size(G);
	X = 0:999;
	Mutations = asLogical(rbinom(1000,1,0.1));
	for (i in X){
		if(Mutations[i]){
			for (g in G){	
				g.addNewMutation(m1,1.0,i);
			}
		}
	}	
	
	catn(size(sim.mutations.position));		
	catn(unique(sim.mutations.position));

}
// output samples of 10 genomes periodically, all fixed mutations at end
1000 late() { p1.outputSample(10); }
2000 late() { p1.outputSample(10); }
2000 late() { sim.outputFixedMutations(); }
