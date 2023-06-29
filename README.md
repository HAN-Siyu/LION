# LION: An Integrated R Package for Effective Prediction of <ins>L</ins>ncRNA- and ncRNA-Prote<ins>I</ins>n Interacti<ins>ON</ins>

Understanding ncRNA-protein interaction is of critical importance to unveil ncRNAs' functions. Now many computational tools have been developed to facilitate the research on ncRNA-protein interaction. Nonetheless, the majority of these tools show unstable results and lack the flexibility required by dataset-specific prediction. Here we propose an integrated package LION which comprises a new method for predicting ncRNA/lncRNA-protein interaction as well as a comprehensive strategy to meet the requirement of customisable prediction. As an integrated tool for predicting ncRNA-protein interaction, LION can be used to build adaptable models for species and tissue-specific prediction and considerably enhance the performance of several widely-used tools. Experimental results also demonstrate our method outperforms its competitors on multiple benchmark datasets. We expect LION will be a powerful and efficient tool for the prediction and analysis of ncRNA- and lncRNA-protein interaction.

## Install LION

**Using devtools**


```
# Enter the following command in R:

if (!library("devtools", logical.return = T)) install.packages("devtools")
devtools::install_github("HAN-Siyu/LION")
```

**Or Download Source Package [Here](https://github.com/HAN-Siyu/LION_Supplementary/raw/master/LION_0.2.9.1.tar.gz) and Install Manually.**

Versions below v0.2.9.1 has a issue in calculating metrics. The issue did not affect the results reported in our paper. We recommend using the latest version. Update details can be found in NEWS.

## Supporting Files

[[PDF Manual](https://github.com/HAN-Siyu/LION_Supplementary/blob/master/LION_0.2.9.1.pdf)]
[[Datasets and Raw Results](https://github.com/HAN-Siyu/LION_Supplementary)]

## Dependencies

Almost all dependencies have been installed when installing LION. However, secondary strucutre features are computed using standalone software, RNAsubopt (from ViennaRNA package) and Predator. You need to download these two programmes if you would like to use method lncPro or extract structural features.

* ViennaRNA Package: https://www.tbi.univie.ac.at/RNA/
* Predator: https://bioweb.pasteur.fr/packages/pack@predator@2.1.2

## Basic Guideline

We expect LION could be a powerful package for predicting RNA-protein interaction in a uniform R environment. The functions of LION can be categorized into several groups to facilitate feature extraction, interaction prediction and model tuning. We here provide a basic summary for LION's function. Detailed examples and parameters explanations can be found in our [manual](https://github.com/HAN-Siyu/LION_Supplementary/blob/master/LION_0.2.8.pdf).

**Functions for feature extraction**
- `computeFreq()`: compute *k*-mer frequencies of RNA/protein sequences. Support three amino acids reprentations, entripy density profile (EDP) computation and data normalization.
- `computeMLC()`: compute the most-like coding region of RNA sequences. Support two strategies: longest open reading frame (ORF) and maximum subarray sum (MSS).
- `computeMotifs()`: compute number of motifs in RNA/protein sequences. User-defined motifs are also supported.
- `computePhysChem()`: compute physicochemical features of RNA/protein sequences. See the [manual](https://github.com/HAN-Siyu/LION_Supplementary/blob/master/LION_0.2.8.pdf) for details.
- `computePhysChem_AAindex()`: compute various physicochemical features of protein sequences using [AAindex](https://www.genome.jp/aaindex/aaindex_help.html).
- `computeStructure()`: computes the secondary structural features of RNA/protein
sequences using [ViennaRNA](https://www.tbi.univie.ac.at/RNA/index.html)/[Predator](https://bioweb.pasteur.fr/packages/pack@predator@2.1.2) packages (the packages are required).

**Functions for feature set construction**
- `featureFreq()`: calculate and construct feature set using *k*-mer frequencies.
- `featureMotifs()`: calculate and construct feature set using motif patterns.
- `featurePhysChem()`: calculate and construct feature set using physicochemical properties.
- `featureStructure()`: calculate and construct feature set using the secondary structural information.

**Functions for random forestion model training**
- `randomForest_CV()`: perform stratified *k*-cross-validation.
- `randomForest_RFE()`: perform stratified feature selection using recursive feature elimination (RFE).
- `randomForest_tune()`: tuning `mtry` of random forest model.
 
**Functions for RNA-protein prediction with different methods**
- `run_LION()`: predict interaction or construct feature set or retrain models using LION method (this work).
- `run_LncADeep()`: predict interaction (retrained random forest model) or construct feature set or retrain models using [LncADeep](https://academic.oup.com/bioinformatics/article/34/22/3825/5021677) method. If you would like to use original deep neural network-based model, please refer to the original [repository](https://github.com/cyang235/LncADeep).  
- `run_lncPro()`: predict interaction (support original algorithm and retrained random forest model) or construct feature set or retrain models using [lncPro](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-651) method. [Original repository](http://cmbi.bjmu.edu.cn/lncpro) is not available when publishing this readme document.
- `run_rpiCOOL()`: predict interaction (retrained random forest model) or construct feature set or retrain models using [rpiCOOL](https://www.sciencedirect.com/science/article/abs/pii/S0022519316300534) method. [Original repository](http://biocool.ir/rpicool.html) is not available when publishing this readme document.
- `run_RPISeq()`: predict interaction (support [web-based original algorithm](http://pridb.gdcb.iastate.edu/RPISeq/) and retrained random forest model) or construct feature set or retrain models using [RPISeq](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-489) method. 
- `run_confidentPrediction()`: perform confident prediction by employing all available methods. Users can further calculate intersection/union or build new models with the output of this function.

**Other Utilities**
- `formatSeq()`: generate sequences pairs for feature extraction or prediction.
- `evaluatePrediction()`: compute metrics, including TP, TN, FP, FN, Sensitivity, Specificity, Accuracy, F1-Score, MCC (Matthews Correlation Coefficient) and Cohen’s Kappa, to evaluate prediction results.
- `runPredator()`: call Predator to process protein sequences ([Predator](https://bioweb.pasteur.fr/packages/pack@predator@2.1.2) is required).
- `runRNAsubopt()`: call RNAsubopt to process protein sequences ([ViennaRNA package](https://www.tbi.univie.ac.at/RNA/index.html) is required).

## Citation

To cite LION in publications, please use:

Siyu Han, Xiao Yang, Hang Sun, Yang Hu, Qi Zhang, Cheng Peng, Wensi Fang, Ying Li. LION: an integrated R package for effective ncRNA-protein interaction prediction. (submitted)


The authors would be glad to hear how LION is used in your study. You are kindly encouraged to notify us (siyu.han@tum.de) about any work you publish!
