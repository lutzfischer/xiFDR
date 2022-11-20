xiFDR
=====

xiFDR is an application for estimating false discovery rates (FDRs) in crosslinking mass spectrometry. It filters crosslinked peptide spectra matches (CSMs) to a list of identifications and derives associated confidence values (False Discovery Rate).

It performs generic FDR calculations for CSMs and resulting peptide pairs, crosslinks and protein pairs.

You can download the latest release of xiFDR from 
https://www.rappsilberlab.org/software/xifdr/ . xiFDR is implemented as a java application and requires java 8 or above to run.

### Background
Correct estimation of false discovery rates in crosslinked peptide identifications presents several quirks that must be taken into account. The 2 main issues are 1) correct estimation of FDR from the target-decoy approach at the level of CSMs and 2) propagation of error from CSMs to peptide pairs, crosslinked residue pairs and protein-protein interactions. xiFDR handles both of these issues allowing for accurate FDR estimation from a search result file.

1. Estimation of FDRs in crosslinking MS via the target-decoy approach

The theoretical basis of adapting the proteomic target-decoy approach for crosslinking MS are covered in [Waltzoheni et al. 2012](https://www.nature.com/articles/nmeth.2103) and [Fisher et al. 2017](https://pubs.acs.org/doi/pdf/10.1021/acs.analchem.6b03745)  and updated for heterobifunctional crosslinkers (crosslinkers with different reactivities at either end) in [Fisher et al. 2018](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0196672). The approach outlined in those papers is implemented in xiFDR. Briefly, decoy protein sequences are added to the database and spectra can then be matched to target linear peptides (T), decoy linear peptides (D) or to crosslinked peptide pairs made up of 2 target sequences (TT), a target-decoy (TD) sequence pair or a decoy-decoy pair (DD).

The formula to estimate FDR for CSMs is then

    FDR=(TD-DD)/TT

xiFDR performs this calculation by finding the score cutoff that corresponds to a desired FDR number at the level of analysis (see next section). 

2. Self and heteromeric crosslinks



3. Error propagation

Very often, the final product of a crosslinking MS experiment is not a list of spectral matches, but rather a list of which residues were found crosslinked to which other residues ("crosslinks" or "links" or "crosslinked residue pairs"). To obtain this list, CSMs must be aggregated into which peptides are paired with each other (multiple spectra can come from the same peptide pair), and those need to be aggregated into crosslinked residues (multiple peptide pairs can come from the same crosslinked residues, because of modifications and miscleavages, for example). Finally, if we are interested in producing a protein-protein interaction (PPI) network, error has to be propagated from residue pairs to PPIs. Correct error propagation from lower to higher levels of result aggregation has big implications for the error rate of the final 
result, as
covered in [Lenz et al. 2021](https://www.nature.com/articles/s41467-021-23666-z) and [Yugandhar et al. 2020](https://www.nature.com/articles/s41592-020-0959-9). Thus, in reporting crosslinking MS data on which residues are crosslinked to each other, FDR filtering should be set at the link level and not at the CSM level. 

xiFDR allows the user to filter for the desired FDR at the level of interpretation of the results. For example, data may be filtered at 5% FDR at the link level, and 10% at the PPI level. Error is propagated by aggregating target and decoy matches from lower levels with a sum of squares approach.


#### Terminology
| Term                        | Description |
|-----------------------------| ----------- |
| CSM                         | Title       |
| PSM                         | Title       |
| peptide pair                | Title       |
| Link/crosslink/residue pair | Title       |
| PPI                         | Title       |
| linear peptide              | Title       |
| Self link                   | Title       |
| Heteromeric link            | Title       |
| Ambiguity                   | Title       |
| Protein Group               | Title       |
| Prefilter                   | Title       |
| Local FDR                   | Title       |
| Posterior Error Probability | Title       |
| Boosting                    | Title       |
| DeltaScore                  | Title       |
| Conservative                | Title       |
| Coverage                    | Title       |




### The interface
The interface provides several tabs. 

The first tab, "input", is used to read in results from search engines, define column names if necessary, and apply prefilters prior to FDR estimation.

The second, "FDR settings", performs the actual FDR calculation at the desired error level.

The "Results" tab shows a summary of the FDR calculation. It also allows the user to save FDR-filtered crosslinking MS results in .csv or [HUPO-compliant mzIdentML 1.2.0 format](https://www.psidev.info/mzidentml#mzid12) for deposition in databases. The output can also be directly uploaded to [xiView.org](https://xiview.org/xiNET_website/index.php) for spectral analysis, network visualization and mapping to structures.

There is then a tab displaying a log of the FDR calculation and an "about" tab with software information.

#### Loading search results from xiSEARCH
xiSEARCH outputs results in a .csv file that contains all crosslinked peptide spectrum matches (CSMs) and scores. It includes target-target matches, target-decoy and decoy-decoy matches.

#### Loading search results from other crosslinking MS search engines

#### Performing FDR filtering

#### FDR suggestions & gotchas

#### Boosting

#### Writing search results

### Custom FDR settings and filters

### Running xiFDR from the command line