## LION: An Integrated R Package for Effective Prediction of <ins>L</ins>ncRNA- and ncRNA-Prote<ins>I</ins>n Interacti<ins>ON</ins>

Understanding ncRNA-protein interaction is of critical importance to unveil ncRNAs' functions. Now many computational tools have been developed to facilitate the research on ncRNA-protein interaction. Nonetheless, the majority of these tools show unstable results and lack the flexibility required by dataset-specific prediction. Here we propose an integrated package LION which comprises a new method for predicting ncRNA/lncRNA-protein interaction as well as a comprehensive strategy to meet the requirement of customisable prediction. As an integrated tool for predicting ncRNA-protein interaction, LION can be used to build adaptable models for species and tissue-specific prediction and considerably enhance the performance of several widely-used tools. Experimental results also demonstrate our method outperforms its competitors on multiple benchmark datasets. We expect LION will be a powerful and efficient tool for the prediction and analysis of ncRNA- and lncRNA-protein interaction.

**Install LION**:

```
if (!library("devtools", logical.return = T)) install.packages("devtools")
devtools::install_github("HAN-Siyu/LION")
```
You can find the relevant data (including datasets and raw results) of this study at https://github.com/HAN-Siyu/LION_Supplementary

**Dependencies**:

Almost all dependencies have been installed when installing LION. However, secondary strucutre features are computed using standalone software, RNAsubopt (from ViennaRNA package) and Predator. You need to download these two programmes if you would like to use method lncPro or extract structural features.

* ViennaRNA Package can be downloaded from: https://www.tbi.univie.ac.at/RNA/
* Predator can be found at: https://bioweb.pasteur.fr/packages/pack@predator@2.1.2
