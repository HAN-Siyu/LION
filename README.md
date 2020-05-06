# ncProR: An Integrated R Package for Effective Prediction of ncRNA- and lncRNA-protein Interaction

Understanding ncRNA-protein interaction is of critical importance to unveil ncRNAs' functions. Now many computational tools have been developed to facilitate the research on ncRNA-protein interaction. Nonetheless, the majority of these tools show unstable results and lack the flexibility required by dataset-specific prediction. Here we propose an integrated package ncProR which comprises a new method for predicting ncRNA/lncRNA-protein interaction as well as a comprehensive strategy to meet the requirement of customisable prediction. As an integrated tool for predicting ncRNA-protein interaction, ncProR can be used to build adaptable models for species and tissue-specific prediction and considerably enhance the performance of several widely-used tools. Experimental results also demonstrate our method outperforms its competitors on multiple benchmark datasets. We expect ncProR will be a powerful and efficient tool for the prediction and analysis of ncRNA- and lncRNA-protein interaction.

*Install ncProR*:

```
if (!library("devtools", logical.return = T)) install.packages("devtools")
devtools::install_github("HAN-Siyu/ncProR")
```
You can find the relevant data (including datasets and raw results) of this study at https://github.com/HAN-Siyu/ncProR_Supplementary