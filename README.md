# Hbond
Python script for identifying presence of hydrogen bonds based on geometrical criteria

## General description

Quantitative analysis of hydrogen bonds is far from trivial. The cause of this is that there is no direct experimental observable for hydrogen bonds.
This has a consequence that numerous ways of defining the hydrogen bond appear in the literature, where each of these definitions emphasizes different aspects of hydrogen bond.
Criteria used to define the presence of hydrogen bond can be divided into energy- and geometry-based criteria. However, this is not the most precise classification, because using
geometric criteria does not mean that energy aspects are excluded. They are to some extent implicitely contained in stereochemistry of the hydrogen bond.
Even though the geometric criteria ofter overestimate the occurence of hydrogen bonds, they are far more computationally efficient and convenient compared to energy-based analysis.
However, finding the most ideal geometric definition of hydrogen bond is not an easy task, and moreover they often contain certain degree of arbitrarity and approximation.

In this script I applied slightly modified geometric definition of that proposed by Raschka _et al_, 2018.<sup>1</sup> This definition of hydrogen bond analyzes local geometry of four atoms, namely donor, hydrogen, acceptor and pre-acceptor, 
which are user-defined according to mdtraj selection syntax, as follows:

```
donor = "resSeq 15 and name NZ"
hydrogen = "resSeq 15 and name HZ1"
acceptor = "resSeq 111 and name OE1"
pre_acceptor = "resSeq 111 and name CD"
```

For hydrogen bond to be identified as "occurred", three conditions are required to be simultaneously satisfied:

* distance between hydrogen and acceptor to be **d<sub>H-A</sub> &le; 2.5 Å**
* angle between donor, hydrogen and acceptor to be in range **120° &le; &theta; <sub>D-H-A</sub> &ge; 180°**
* angle between hydrogen, acceptor and pre-acceptor atom to be in range **90° &le; &theta; <sub>H-A-PA</sub> &ge; 180°**

These three geometric criteria are atom-specific, and in principle should be tailored for a given combination of atoms. d<sub>H-A</sub> can be estimated from the first minimum of the radial distribution function of a 
given pair of atoms, while angles  &theta; <sub>D-H-A</sub> and &theta; <sub>H-A-PA</sub> can be obtained from X-ray measurements.

For each frame in trajectory, an observable _h_ takes value 1 if all three criteria are met and hydrogren bond occurs. Otherwise, it takes value zero.
Occurence of hydrogen bond, as an array of ones and zeros, is saved as .dat file. 

**Outlook**

This binary time trace of hydrogen bond occurence can be subsequently analyzed. For instance, one can calculate so-called continuous hydrogen bond autocorrelation function, which can be integrated to find the hydrogen bond lifetime.<sup>2</sup>
Since presence of hydrogen bond can be a property of enzymatically active or inactive species, another possibility is to therefore use hydrogen bond occurence for clustering/slassification of structural models.

## Input file requirements

* all-atom MD trajectory in any of the mdtraj compatible formats (dcd, nc, xtc..)
* topology as .pdb file


## Dependencies
_Hbond_analysis.py_ is a python script built on Python 3.8.8. Script was tested under the following configuration:

* Windows 10
* Python 3.8.8
* mdtraj 1.9.4
* numpy 1.23.0


## References

1. Raschka, S., Wolf, A.J., Bemister-Buffington, J. et al.
Protein–ligand interfaces are polarized: discovery of a strong trend for intermolecular hydrogen bonds to favor donors on the protein side with implications for predicting and designing ligand complexes. J Comput Aided Mol Des 32, 511–528 (2018)

2. Ippolito JA, Alexander RS, Christianson DW (1990) Hydrogen bond stereochemistry in protein structure and function. J Mol Biol 215:457–471

3. Richard J. Gowers and Paola Carbone, A multiscale approach to model hydrogen bonding: The case of polyamide The Journal of Chemical Physics, 142, 224907 (2015)



## Authors

* Milana Popara
