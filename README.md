xiFDR
=====

xiFDR is an application for estimating false discovery rates (FDRs) in crosslinking mass spectrometry. It filters crosslinked peptide spectra matches (CSMs) to a list of identifications and derives associated confidence values.

It performs a generic FDR calculations for CSMs and resulting peptide pairs, crosslinks and protein pairs. It complies with the standards for data reporting set by the HUPO proteomics standards initiative and can output results in .mzIdentML 1.2.0 format for deposition in databases. It is search engine-agnostic and can therefore perform FDR filtering on results from [xiSEARCH](https://github.com/Rappsilber-Laboratory/xisearch) but also other crosslinking MS search engines. The output can then be directly uploaded to [xiView.org](https://xiview.org/xiNET_website/index.php) for spectral analysis, network visualization and mapping to structures.

You can download the latest release of xiFDR from https://www.rappsilberlab.org/software/xifdr/ . xiFDR is implemented as a java application and requires java 8 or above to run.

For questions regarding usage of xiFDR, please open a [discussion](https://github.com/Rappsilber-Laboratory/XiSearch/discussions).

### Background
Correct estimation of false discovery rates in crosslinked peptide identifications presents several quirks that must be taken into account. The 2 main issues are 1) correct estimation of FDR from the target-decoy approach at the level of CSMs 2) Correct handling of the random space for self and heteromeric crosslinks and 3) propagation of error from CSMs to peptide pairs, crosslinked residue pairs and protein-protein interactions. xiFDR handles both of these issues allowing for accurate FDR estimation from a search result file.

1. Estimation of FDRs in crosslinking MS via the target-decoy approach

The theoretical basis of adapting the proteomic target-decoy approach for crosslinking MS are covered in [Waltzoheni et al. 2012](https://www.nature.com/articles/nmeth.2103) and [Fisher et al. 2017](https://pubs.acs.org/doi/pdf/10.1021/acs.analchem.6b03745) and updated for heterobifunctional crosslinkers (crosslinkers with different reactivities at either end) in [Fisher et al. 2018](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0196672). The approach outlined in those papers is implemented in xiFDR. Briefly, decoy protein sequences are added to the database and spectra can then be matched to target linear peptides (T), decoy linear peptides (D) or to crosslinked peptide pairs made up of 2 target sequences (TT), a target-decoy (TD) sequence pair or a decoy-decoy pair (DD).

The formula to estimate FDR (provided the chance of matching decoys accurately models the chance of random matching) for CSMs is then

    FDR=(TD-DD)/TT

xiFDR performs this calculation by finding the score cutoff that corresponds to a desired FDR number at the level of analysis (see next section). 

2. Self and heteromeric crosslinks

The chance of random matching a crosslinked peptide pair spectrum within a protein sequence (self crosslink) is different from that of matching a crosslink between 2 different proteins (heteromeric crosslink) as shown in [Lenz et al. 2021](https://www.nature.com/articles/s41467-021-23666-z) Fig. 1. Moreover, spectra of heteromeric crosslinks tend to have lower signal-to-noise than those of self crosslinks. For these reasons, FDR calculation on self and heteromeric crosslinks should be performed separately. xiFDR automatically splits FDR calculations between self and heteromeric crosslinks reporting them separately in the "results" page. 

3. Error propagation



Very often, the final product of a crosslinking MS experiment is not a list of spectral matches, but rather a list of which residues were found crosslinked to which other residues ("crosslinks" or "links" or "crosslinked residue pairs"). To obtain this list, CSMs must be aggregated into which peptides are paired with each other (multiple spectra can come from the same peptide pair), and those need to be aggregated into crosslinked residues (multiple peptide pairs can come from the same crosslinked residues, because of modifications and miscleavages, for example). Finally, if we are interested in producing a protein-protein interaction (PPI) network, error has to be propagated from residue pairs to PPIs. Correct error propagation from lower to higher levels of result aggregation has big implications for the error rate of the final 
result, as
covered in [Lenz et al. 2021](https://www.nature.com/articles/s41467-021-23666-z) and [Yugandhar et al. 2020](https://www.nature.com/articles/s41592-020-0959-9). Thus, in reporting crosslinking MS data on which residues are crosslinked to each other, FDR filtering should be set at the link level and not at the CSM level.

![ErrorPropagation](https://user-images.githubusercontent.com/44289027/203276815-046f76bf-2c1d-49a0-ac7a-479c73b32a3a.PNG)

(target-target matches in grey, matches involving decoys in red. Error control should account for aggregation of matches from lower levels to higher levels.)

xiFDR allows the user to filter for the desired FDR at the level or levels of interpretation of the results. For example, data may be filtered at 5% FDR at the link level, and 10% at the PPI level. Error is propagated by aggregating target and decoy matches from lower levels with a sum of squares approach.


#### Terminology
| Term                                     | Description                                                                                                                                     |
|------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| CSM                                      | crosslink spectra match                                                                                                                         |
| PSM                                      | peptide (or crosslink) spectra match                                                                                                            |
| peptide pair                             | 2 peptides crosslinked to each other, no matter the crosslink site. Multiple CSMs can support the same peptide pair                             |
| Link/crosslink/residue pair              | 2 residues crosslinked to each other. Multiple peptide pairs can support the same crosslink                                                     |
| Protein pair/PPI                         | 2 proteins found crosslinked to each other, no matter the number/identity of the residue pairs. Multiple residue pairs can support the same PPI |
| linear peptide                           | a non-crosslinked peptide                                                                                                                       |
| Self link                                | a link within a protein sequence (includes links between multiple copies of the same protein)                                                   |
| Heteromeric link                         | a link between 2 different protein sequences                                                                                                    |
| Ambiguity                                | when peptide pairs or residue pairs can be mapped to multiple proteins because these share the same sequence                                    |
| Protein Group                            | the group of proteins sharing the sequence making up the ambiguous group                                                                        |
| Prefilter                                | a score filter applied prior to FDR calculation to both target and decoy matches                                                                |
| Local FDR/poterior error probability/PEP | the estimation of FDR of a particular match based on a windowed approach or a fitting of the score distributions.                               |
| Boosting                                 | a grid search approach optimising settings to maximise the number of matches passing a given FDR threshold                                      |
| DeltaScore                               | The score of the best explanation of a CSM/residue pair etc. divided by the second best explanation                                             |
| Conservative                             | An explanation of a spectral feature where non-lossy matches are weighted more heavily than lossy ones.                                         |
| Coverage                                 | number of fragments matched / max number of theoretical fragments                                                                               |

## Calculating FDR for crosslinking MS data with xiFDR

### The interface
The interface should be loaded via the executable files "startWindows.bat"/"startUnix.sh"/"startMacOS.command" depending on whether the program is run in Windows, Linux or Mac. It is not advisable to run the .jar file directly as this may not detect the correct encoding, which will impact writing the .mzIdentML file.

The interface provides several tabs. 

The first tab, "input", is used to read in results from search engines, define column names if necessary, and apply prefilters prior to FDR estimation.

The second, "FDR settings", performs the actual FDR calculation at the desired error level.

The "Results" tab shows a summary of the FDR calculation. It also allows the user to save FDR-filtered crosslinking MS results in .csv or [HUPO-compliant mzIdentML 1.2.0 format](https://www.psidev.info/mzidentml#mzid12) for deposition in databases.

There is then a tab displaying a log of the FDR calculation and an "about" tab with software information.

### Loading crosslinking MS search results

#### Loading search results from xiSEARCH
xiSEARCH outputs results in a .csv file that contains all crosslinked peptide spectrum matches (CSMs) and scores. It includes target-target matches, target-decoy and decoy-decoy matches.

Select the output .csv file. The .fasta database and the .config file used in the search should also be uploaded in the respective boxes if one plans to generate an .mzIdentML result file.

At this stage, any prefilters on spectral quality (see below) should be set by ticking the "filter" box.

Press "read" to initiate reading of result. Monitor  the bottom left of the window for the message "finished reading file" and the bottom right for memory usage. 

If xiFDR is very slow or crashing, restart the program by increasing the allocated memory by editing the java -Xmx option in the startup .bat/.sh/.command file. The default is -Xmx3G, providing 3 Gb of RAM. Often, this will not be enough for searches with dozens of runs and thousands of matches. Large searches may require tens of Gb to read in. In this case, prefilters (see below) are a useful way to reduce memory requirements if needed.


#### prefilters
xiSEARCH provides many features of CSMs that may be used to prefilter the results prior to FDR estimation. Doing so equally on targets and decoys prior to FDR estimation retains the accuracy of FDR, while doing score filtering or other filters post FDR estimation generates results with unknown error rates.

The prefilters may be toggled in the "input" tab by clicking the "filter" option. These are generally used after a first look at results without prefilters if spectra of low quality are still passing the FDR. Some of the commonly used prefilters (for FDR calculation performed on xiSEARCH results, especially on large scale searches) are

| Filter name                                         | Description                                                                                                                                                    | Commonly set to               |
|-----------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------|
| peptide 1/2 unique matched                          | minimum observed fragments for peptide 1 or 2                                                                                                                  | >3 (4-5 to be more stringent) |
| fragment unique crosslinked matched conservative    | minimum fragments across both peptides containing the crosslink site                                                                                           | >0                            |
| peptide 1/2 unique crosslinked matched conservative | impose to have only spectra where both peptides have a fragment of the other                                                                                   | >0                            |
| delta                                               | only spectra where the best match is x times better than the second best by score                                                                              | >1.2                          |
| peptide1/2 CCPepFragmentDoubletCount                | For cleavable crosslinkers, only consider spectra where at least X doublets are found on peptide 1/2                                                           | >0                            |
| fragment CCPepDoubletCount                          | For cleavable crosslinkers, only consider spectra where at least X doublets are found on either peptide. usually set to greater than 0 in large scale searches | >0                            |
Notice that these are not meant to be used blindly all at once!

MS-cleavable crosslinkers present several advantages. Their signature crosslinker stubs and peptide doublets help to provide extra confidence, that the peptide masses are correct - and in extend increases the chance that the the peptides themself are correctly identified. xiFDR can make the most out of these features by prefiltering spectra on a minimum of crosslinker stubs observed, and then boosting on stubs and doublets.


#### Loading search results from other crosslinking MS search engines
xiFDR can directly read search result files in te standard .mzIdentML format in the mzIdentML tab. 

If the results are not in .mzIdentML format, search results should be read in via the csv tab. To do so, a minimal csv has to be generated containing the following information:


| column name | Description |
|-------------| ----------- |
| XXX         | Title       |
| XXX         | Title       |
| XXX         | Title       |
| XXX         | Title       |
| XXX         | Title       |

Below are a few tips for specific search engines:
- The decoys have to be reported in the search results and they should be then mapped to target proteins
- XXX


Column names may be remapped in the bottom half of the csv interface by XXXX

#### Loading results from ProteomeXchange repositories

#### other settings in input tab

| setting      | Description                                                                                                             | When to use                                                                                                                                                                                                                   |
|--------------|-------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Flag modifications | splits FDR calculation for modified and unmodified peptides                                                             | if spectral quality of modified and unmodified peptides differs and there are sufficient numbers of matches in both modified and unmodified to perform FDR calculation                                                        |
| Normalize    | Normalize score distributions according to FDR or the decoy distribution. Suggested settings "FDR based" or "AutoScore" | generate results to be merged with other searches with different settings. For example, normalize scores to later combine FDR results from searches with HCD and ETD fragmentation, that yield different score distributions. |
| Forward Columns    | Text box for a regular expression to take columns from the input .csv file and include it in the xiFDR output files | include for example retention time in the result file to perform additional rescoring based on machine learning, for example with [xiRT](https://www.nature.com/articles/s41467-021-23441-0)|



### Performing the FDR calculation

In the "FDR settings", one can perform the actual FDR filtering. 

By default, the view is set to "reduced FDR", which shows just the basic settings. The cutoff is set at 5% at the residue pair level, meaning the  error will be propagated so that 5% of the residue pairs correspond to a wrong/random match. The "boosting" features is enabled (see below for more details). These are perfectly acceptable FDR filtering settings for experiments aimed at characterising the crosslinks in a purified protein complex and should give a good idea of the number of crosslinks detectable with reasonable certainty in the sample. In analyses of cellular fractions or searches with hundreds of proteins in the database, it is advisable to also include an FDR cutoff at the "Protein Pairs" level. Similarly, in analyses devoted to method development on the quality of spectra, a filter at the "PSM" is advised.

If these settings are satisfactory, press "calculate". 

Otherwise, more advanced options are available if one ticks the "complete FDR" setting. Further ticking the box "more" will present the full set of options for FDR filtering in xiFDR. Hovering the mouse over checkboxes or up/down arrows provides a brief explanation of the setting.

All FDR settings:

| Setting                                                     | Description                                                                                                                                                                                                                                        | when to use                                                                                                           |
|-------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------|
| max FDR                                                     | the FDR setting at that particular level (PSM, peptide pair etc.)                                                                                                                                                                                  |                                                                                                                       |
| Local FDR                                                   | whether to compute a per-PSM/per-peptide pair etc. score based on a windowed FDR approach. Black box=compute local FDR, ticked box=calculate and filter by local FDR, white box=do not calculate local FDR                                         | to analyse distributions of scores and FDR or to generate a per-crosslink error estimate e.g. for structure modeling. |
| min Pep. length                                             | minimum length of peptide 1 and peptide 2                                                                                                                                                                                                          | increase to 6 if noise matching of short sequences                                                                    |
| min supporting peptide pairs Protein Group/Residue Pair/PPI | minimum number of peptide pairs necessary to aggregate to a given level                                                                                                                                                                            | define to have more robust network especially at the PPI level.                                                       |
| min TD chance                                               | require enough target matches so that these number of decoy matches would still pass the target FDR. In other words, a minimum number for reliable FDR calculation                                                                                 | Use for more accurate FDR calculation. By default, this number is ignored because "ignore validity checks" is ticked  |
| Ignore validity checks                                      | Allow FDR calculation even if number of decoys is not sufficient for sufficiently accurate estimation of FDR                                                                                                                                       | Untick and use in combination with min TD chance to only get results with accurate FDR calculation                    |
| Ambiguity                                                   | increase number to only allow matches with at most N ambiguity. E.g. if "ambiguity" is set to 1 in the residue pairs, residue pairs that could come from more than 1 place in the sequence or more than 2 unique proteins are discarded            | Set to a number if not interested in considering ambiguous crosslinks at all |
| Unique PSMs                                                 | Aggregate specral matches so that only the top scoring spectrum is considered per spectral match                                                                                                                                                   | on by default. Untick if interested in an analysis of spectral matches                                                                       |
| Group by PSMs                                               | split FDR calculations into groups of protein pairs having different numbers of matches. In theory splittng FDR for highly abundant proteins and low-abundant                                                                                      | off by default. Experimental                                                                                                                 |
| No consecutive                                              | remove matches from consecutive peptides. 2 consecutive (in sequence) linear peptides of which one is crosslinker modified can have the same mass and even fragmentation pattern of a crosslinked peptide pair and tend to be not very informative | Consider switching on in searches were protein interactions are the focus, or searches where short range sequence contacts are not the focus |
| Unique PSMs                                                 | Aggregate specral matches so that only the top scoring spectrum is considered per spectral match                                                                                                                                                   | on by default. Untick if interested in an analysis of spectral matches                                                                       |

Settings in the "more" panel:

These are minimum filters that essentially act as prefilters prior to FDR calculation. The difference between setting them here versus in the "input" tab is that these parameters can then be part of the boosting grid search (see below).

| Setting | Description| when to use                                                    |
|---------|------------|----------------------------------------------------------------|
| min. peptide fragments| minimum observed fragments per peptide | 0 by default, increase number in SDA searches                  |
| min. coverage| number between 0 and 1. Minimum fraction of theoretical fragments required | included in boosting                                           |
| min. coverage| number between 0 and 1. Minimum fraction of theoretical fragments required | included in boosting                                           |
| min. peptide stubs| minimum number of fragment stubs observed | important for MS-cleavable crosslinkers, included in boosting  |
|min. peptide doublets| minimum number of  peptide doublets observed | important for DSSO searches, included in boosting              |
| boost separately|perform independent optimisation of each parameter| on by default                                                  |

The settings in "define groups" are currently very beta and unsupported. They allow for splitting the dataset further down into custom groups for FDR calculation.

#### Boosting
xiFDR's boosting feature performs a grid search to optimise FDR settings at lower levels to reach the maximum number of matches passing validation at the desired FDR level. For example, enabling boosting for a residue pair-level FDR of 5% will tweak PSM, peptide pair and protein group level FDR cutoffs to maximise the number of residue pairs passing the 5% FDR threshold.

Boosting is performed with a grid search of parameters as described in [Fisher et al. 2017](https://pubs.acs.org/doi/pdf/10.1021/acs.analchem.6b03745). 

The user may control which parameters are part of boosting by changing the selection in the "boosting includes" button. 

The "steps" controls how many steps of the grid search per parameter. The "between" box ensures that boosting is performed to maximise the number of heteromeric residue pairs/PPIs etc. passing FDR rather than the overall number. This is recommended for searches where the goal is to produce a protein-protein interaction network and where large numbers of heteromeric crosslinks are available.

We recommend leaving boosting on and selecting "between" if desired. For experiments with MS-cleavable crosslinkers, we suggest boosting on minimum peptide stubs and minimum peptide doublets.

FDR calculations with boosting enabled can take some minutes to conclude.

### FDR suggestions & common mistakes

The extreme flexibility of xiFDR requires some careful use. Here are some of our suggestions and ways to spot that the FDR calculation may not be giving sensible results. When in doubt, keep it simple! 

#### Tips & tricks
- For searches with cleavable crosslinkers such as DSSO, set prefilters (see below) on the number of peptide doublets and boost on that parameter
- For searches with unspecific crosslinkers like SDA, include prefilters (or boost) on minimum number of fragments. Include also a prefilter on "fragment unique crosslinked matched conservative" and exclude consecutives. This ensures that only spectral matches with at least one fragment proving a crosslinked peptide pair are considered.

#### Watch out for:
- The program warns if there aren't sufficient targets to estimate FDR accurately. If this warning refers to the level of analysis you are interested in, there is likely almost nothing in the data. Unless prefilters were set way too stringent.
- The FDR calculation rests on the assumption that target-decoy pairs have twice the chance of random matching than decoy-decoy pairs (TD+DT vs DD). If in your results you have more decoy-decoy (DD) than target-decoy (TD), your FDR evaluates to a negative number and becomes meaningless. This may be a sign that there are no crosslinks in the datasets, or that prefilters are not working as intended.


### Results summary
After the calculation is complete (check the log tab, the lower bottom left of the window or the "calculate" button becoming clickable again), the results tab presents a summary of all matches passing the FDR threshold with the number of decoys and target-decoys in each category. Remember the ratio of targets and decoys passing determines the accuracy of the FDR estimation.

### Writing out search results

Results may be written out in xiFDR .csv format (csv tab) or in mzIdentML1.2.0 format. We advise to do both every time, as the summary file is useful to keep track of what was done, and the mzIdentML for later deposition or upload to [xiView.org](https://www.xiview.org).

#### csv output 
will generate several files:

- Summary file: a file containing a summary of the FDR calculation, results and all settings
- CSM file: a file containing all crosslink spectra matches passing the threshold. This file can be uploaded to xiview.org for visualization
- peptide pairs file: a file containing all peptide pairs passing the threshold
- protein groups file: all protein groups passing the threshold (i.e. including ambiguity when a peptide cannot be assigned to an individual protein)
- Links file: all the crosslinked residue pairs passing the threshold
- PPI file: all the protein-protein interactions passing the thresold.

Each file contains information about the FDR, local FDR (if enabled), posterior error probability, target/decoy nature of each match, as well as the peakfile of origin.

Beware that if multiple FDR thresholds are set, *only* the matches passing the highest level of aggregation are included in the results. In other words, the csv files for a search with both a 5% residue pair FDR and 5% protein pair FDR will only include residue pairs with a 5% or better FDR that also lead to a protein pair with a 5% or better FDR.

#### mzIdentmL output
will generate a single file .mzIdentML compliant with standards. The file can be deposited in ProteomeXChange repositories or uploaded to xiview.org for visualization. It contains information about the search results, peaks and validation.



### Running xiFDR from the command line

xiFDR may be run in the command line for automated processing or cluster jobs.

The full range of options and their description is available with

    java --jar xiFDR.jar --help

Allocate memory with the -XmX flag.

Prefilters are not supported in the command line version of xiFDR. However, prefiltering the .csv input file as desired may be done in python/pandas, R, or any other tool prior to loading the input file into xiFDR.