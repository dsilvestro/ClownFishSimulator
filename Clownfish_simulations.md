## Clownfish simulations

### Forward simulations with
- 20 colonies of 5 individuals each (female, male, juveniles)
- 1000 nDNA loci for each individual (species A: 0/1 alleles)
- 1 mtDNA locus (species A: 0 or 1)
- recombination frequency at reproduction is set at 3 %
- mortality is set to 1 % of the individuals at each generation. Individuals die with probability varying linearly depending on hierarchical placement (female 3 times more likely to die than male, and 10 times more likely to die than youngest)

### Hybridization events
- at reproduction each colony produces 25 larvae all of which form the larvae pool.
- a fraction (`fB`) of larvae is replaced with larvae of species B (nDNA: 2/3 alleles and mtDNA: 2 or 3 allele).
- larvae have a mortality probability (`delta`) which is a function of their nDNA admixture (`admix` ranging from 0 to 0.5) so that: `delta = (admix*2)*m + min_d`, where `min_d=0.001` is a minimum mortality rate for non-admixed individuals and `m` determines how mortality changes for admixed individuals. if `m + 0.001 = 1` all larvae with maximum admixture (`admix = 0.5`) die.
- dead individuals are replaced by lower-ranking individuals from the same colony and gaps at lowest rank are filled by randomly drawing an individual from larvae pool.

### Parameterization of the simulations
In each simulations the following parameters are sampled:

- number of generations during which species B larvae are added to the pool (hybridization phase): `U[0, 1500]`
- number of generations following the hybridization phase: `U[0, 1500]`
- maximum mortality of hybrid larvae: `m ~ U[0, 1]`
- fraction of species B added to the larvae pool during the hybridization phase: `fB ~ U[0, 1]`

### ABC parameter estimation
Simulations are accepted only if they fulfill the following criteria at the end of the simulation:

- there is at least 1 individual across all colonies with mtDNA from species B
- the median level of admixture among individuals with mtDNA from species B is < 0.2 (i.e. >80% nDNA loci from species A)
- the maximum proportion of nDNA loci from species B across all individuals does not exceed 80%. 