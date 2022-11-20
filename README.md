xiFDR
=====

xiFDR is an application for estimating false discovery rates (FDRs) in crosslinking mass spectrometry. It filters crosslinked peptide spectra matches (CSMs) to a list of identifications and derives associated confidence values (False Discovery Rate).

It performs generic FDR calculations for CSMs and resulting peptide pairs, crosslinks and protein pairs.

You can download the latest release of xiFDR from 
https://www.rappsilberlab.org/software/xifdr/

### Background
Correct estimation of false discovery rates in crosslinked peptide identifications presents several quirks that must be taken into account. The 2 main issues are 1) correct estimation of FDR from the target-decoy approach at the level of CSMs and 2) propagation of error from CSMs to peptide pairs, crosslinked residue pairs and protein-protein interactions. xiFDR handles both of these issues allowing for accurate FDR estimation from a search result file.

1. Estimation of FDRs in crosslinking MS via the target-decoy approach


2. Error propagation

#### Terminology

### The interface
The interface provides several tabs. 

The first tab, "input", is used to read in results from search engines, define column names if necessary, and apply prefilters prior to FDR estimation.

The second, "FDR settings", performs the actual FDR calculation at the desired error level.

The "Results" tab shows a summary of the FDR calculation. It also allows the user to save FDR-filtered crosslinking MS results in .csv or [HUPO-compliant mzIdentML 1.2.0 format](https://www.psidev.info/mzidentml#mzid12) for deposition in databases. The output can also be directly uploaded to [xiView.org](https://xiview.org/xiNET_website/index.php) for spectral analysis, network visualization and mapping to structures.

There is then a tab displaying a log of the FDR calculation and an "about" tab with software information.

#### Loading search results from xiSEARCH

#### Loading search results from other crosslinking MS search engines

#### Performing FDR filtering

#### FDR suggestions

#### Boosting

#### Writing search results

### Custom FDR settings and filters

