# Notes:

Three main things to figure out:

- ways to make the model more realistic
- which parameters to vary/compare to each other
- how to quantitatively compare situations to each other


# Model changes/improvements

1. effect size distribution (not plus/minus 1 but gamma distributed)
2. arrangements of loci along the chromosome ("genomic architecture")


# What to measure

## Clustering of loci along the genome:

How clustered are loci underlying the divergent phenotype clustered together?
We generally want to see if "important" loci tend to be near each other or far from each other,
and possibly to compare:

- freshwater versus saltwater alleles
- large versus small effects
- high and low frequency (fixed, even?)

One way to see if there is clustering would be to compare to a Poisson process
(search "test for uniformity statistics stack exchange Poisson").
Another way to see how much mutations of type A tend to be near ones of type B
(maybe with A=B)
would be: for each mutation of type A, find how many others of type B are within distance x,
as a function of x, 
and divide this by the *expected* number if those of type B were distributed uniformly.
This gives a curve as a function of x, and we can plot these.


## Sharing of adaptation between lakes

How much is the genetic basis of adaptation shared between lakes?
Are haplotypes shared?  Alleles?
Here are some measures:

1. Variance of allele frequencies across lakes.
2. Phenotypic variance of "hybrid" offspring of two parents, one from each of two lakes.
    (related to the Wright-Castle estimator for number of loci)
3. Mean percentage of shared mutations between two fish, as a function of distance apart.



