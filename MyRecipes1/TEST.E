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
	sim.addSubpop("p1", 500);
	
}

2000 late() {
	g = sample(sim.subpopulations.genomes, 10, T);
	m = sortBy(unique(g.mutations), "position"); // print the number of segregating sites
	cat("\n\nsegsites: " + size(m) + "\n");
	positions = format("%.6f", m.position / sim.chromosome.lastPosition);
	cat("positions: " + paste(positions, " ") + "\n");
	for (genome in g) {
		hasMuts = (match(m, genome.mutations) >= 0);
		cat(paste(asInteger(hasMuts), "") + "\n"); 
	}
}
