[![R-CMD-check](https://github.com/PaulodSousa/Haplo-Diploid_Equilibrium/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/PaulodSousa/Haplo-Diploid_Equilibrium/actions/workflows/R-CMD-check.yaml)

# Haplo-Diploid Equilibrium
This repository contains:
1. R scripts to estimate pop. gen. summary statistics from haplo-diploid systems (Script_folder)
2. Fake data (.vcf and .txt) to test scripts (Fake_data)
3. A theoretical explanation of the model (see below)

## List of scripts
* `R/HaploDip_Het_base.R`: Base script where the bases of the model are written
* `R/HaploDip_GenoFreq.R`: Functions to calculate genotype frequencies and F<sub>IS</sub> (not agnostic to ploidy)
* `R/Fst_HapDip.R`: Functions to calculate F<sub>ST</sub> from allele frequency data (agnostic to ploidy)
* `R/Haplo-Dip_Sex_Ref_Allele.R`: Functions to calculate reference allele frequencies (0) in each sex across populations (agnostic to ploidy)
* `R/HaploDip_Neis_H.R`: Functions to calculate Nei's H, a geentic diversity metric (agnostic to ploidy)
* `unfinished/Watternson_Theta.R` > script to calculate Watterson's Theta (agnostic to ploidy) **unfinished**


## Theoretical explanation of the model
In a haplo-diploid system where one sex is haploid and another diploid (e.g. hymentoptera insects) and assuming
1. equal sex ratio
2. random mating
3. infinite population size (no drift)
4. equal probability of reproduction (no selection),
5. alleles don`t change from individual to gene pool (no mutation)
6. no addition of gametes from other gene pools (no gene flow)
7. equal allele frequencies across sexes
8. no overlapping generations

Let p be the frequency of the allele "A" and q be the frequency of the allele "a" of a bi-allelic marker so that p + q = 1, then: \
f(AA) = p² x 1/2 = p²/2 \
f(aa) = q² x 1/2 = q²/2 \
f(Aa) = 2pq x 1/2 = pq \
f(A) = p x 1/2 = p/2 \
f(a) = q x 1/2 = q/2 

f(AA) + f(Aa) + f(aa) + f(A) + f(a) = p²/2 + pq + q²/2 + p/2 + q/2 = \
= 1/2 x (p² + 2pq + q² + p + q) = \
= 1/2 x (p² + 2pq + q²) + 1/2 x (p + q) = \
= 1/2 x (p + q)² + 1/2 x (1) = \
= 1/2 x (1)² + 1/2 = \
= 1/2 + 1/2 = \
= 1

Let N(X) be the number of individuals with genotype X, and N.dip and N.haplo being the number of diploid and haploid individuals respectively, then: \
p = ( 2x N(AA) + N(Aa) + N(A) ) / (2x N.dip + N.haplo) \
q = ( 2x N(aa) + N(Aa) + N(a) ) / (2x N.dip + N.haplo

In order to calculate genotype frequencies on next generation: \
if offspring is diploid, then: \
| Mating pair | Offspring genotype distribution | Mating pair frequency |
|-------------|---------------------------------|-----------------------|
| AA x A      |       AA (100%)                 |  2x ( f(AA) x f(A) )  |
| AA x a      |       Aa (100%)                 |  2x ( f(AA) x f(a) )  |
| Aa x A      |       Aa (50%) + AA (50%)       |  2x ( f(Aa) x f(A) ) |
| Aa x a      |       Aa (50%) + aa (50%)       |  2x ( f(Aa) x f(a) ) |
| aa x A      |       Aa (100%)                 |  2x ( f(aa) x f(A) ) |
| aa x a      |       aa (100%)                 |  2x ( f(aa) x f(a) ) |

if offspring is haploid, then: \
| Parent genotype | Offspring genotype distribution | Parent frequency |
|-----------------|---------------------------------|------------------|
| AA              |       A (100%)                 |       f(AA)      |
| Aa              |       A (50%) + a (50%)         |       f(Aa)      |
| aa              |       a (100%)                 |       f(aa)      |

The frequency of each genotype of next generation are: \
f(AA)t+1 = 2x (f(AA) x f(A)) + 2x (f(Aa) x f(A))/2 = 2x (p²/2 x p/2) + (pq x p/2) = 2x p³/4 + p²q/2 = p³/2 + p²q/2 = p²/2 x (p + q) = p²/2 x 1 = p²/2

f(aa)t+1 = 2x (f(aa) x f(a)) + 2x (f(Aa) x f(a))/2 = 2x (q²/2 x q/2) + (pq x q/2) = 2x q³/4 + q²p/2 = q³/2 + q²p/2 = q²/2 x (q + p) = q²/2 x 1 = q²/2

f(Aa)t+1 = 2x (f(AA) x f(a)) + 2x (f(Aa) x f(a))/2 + 2x (f(Aa) x f(A))/2 + 2x(f(aa) x f(a)) = 2x(p²/2 x q/2) + (pq x q/2) + (pq x p/2) + 2x(q²/2 + p/2) =
= p²q/2 + pq²/2 + p²q/2 + pq²/2 = 1/2 x (p²q + pq² + p²q + pq²) = 1/2 x (2p²q + 2pq²) = 2pq/2 x (p + q) = pq x (1) = pq

f(A)t+1 = f(AA) + f(Aa)/2 = p²/2 + pq/2 = p/2 x (p + q) = p/2 x (1) = p/2

f(a)t+1 = f(aa) + f(Aa)/2 = q²/2 + pq/2 = q/2 x (p + q) = q/2 x (1) = q/2

| Genotype | Frequency at t0 | Frequency at t+1 |
|----------|-----------------|------------------|
| AA       |      p²/2       |        p²/2      |           
| Aa       |      pq         |        pq        |
| aa       |      q²/2       |        q²/2      |
| A        |      p/2        |        p/2       |
| a        |      q/2        |        q/2       |

Since, the genotype frequencies in t+1 are the same as in t0, no evolution as occured.


### if inbreeding is present:
Let f be the probability of a individual's allele be identical by descent with its other allele and (1-f) the probability that is not. Let XX' be the genotype frequency of XX under inbreeding

AA' = fp + (1-f) x p²/2 = fp + p²/2 -fp²/2 = p²/2 + fp(1-p/2)
aa' = fq + (1-f) x q²/2 = fq + p²/2 -fq²/2 = q²/2 + fq(1-q/2)
A' = p/2 and a' = q/2
Aa' = 1 - (AA'+ aa' + A' + a') = 1 - [p²/2 + fp(1-p/2) + q²/2 + fq(1-q/2) + p/2 + q/2] = 1 - [p²/2 + q²/2 + p/2 + q/2 + fp + fq - fp²/2 - fq²/2] = 1 - [(1 - pq) + f(-p²/2 -q²/2 + p + q)] = 1 - [1- pq + f(1 - p²/2 - q²/2)] = pq - f(1 - p²/2 - q²/2)


| Genotype | Frequency with inbreeding |
|----------|-----------------|
| AA       |      p²/2 + fp(1-p/2)       |        
| Aa       |      pq - f(1- p²/2 -q²/2)        |
| aa       |      q²/2 + fq(1-p/2)       |
| A        |      p/2        |
| a        |      q/2        |


| Genotype | f = 0 | f = 1 |
|----------|-----------------|-|
| AA       |      p²/2 + 0       | p |           
| Aa       |      pq - 0        | -p/2 -q/2 |
| aa       |      q²/2 + 0       | q |
| A        |      p/2        | p/2 |
| a        |      q/2        | q/2 |










 
