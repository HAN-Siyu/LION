# Next Release

* Users can use `run_multiple()` to get prediction results of multiple methods in one run.
* Models will be rebuilt using updated data.

# ncProR 1.1.3

* Add `run_LncADeep()` for prediction using LncADeep's feature set.
* Add `computeMLC()` to compute the most-like CDS (MLC) region.
* Add parameter `EDP` which supports applying entropy density profile to sequence intrinsic features.
* Update reference manual.
* Bugs fixed.

# ncProR 1.1.2

* Matthews correlation coefficient (MCC) is included into evaluation metrics.
* TP, TN, FP and FN are reported in evaluation results.
* Add `evaluatePrediction()` for prediction evaluation.
* Bugs fixed.

# ncProR 1.1.1

* Add argument verification.
* Update reference manual.
* Bugs fixed.

# ncProR 1.1.0

* `"retrain"` and `"feature"` modes are added to function `run_ncProR()`, `run_RPISeq()`, `run_lncPro()` and `run_rpiCOOL()`. Now users can use ready-to-use functions to extract features or retrain the models of ncProR, RPISeq, lncPro and rpiCOOL.
* External parallel cores can be passed to functions.
* Internal changes to improve efficiency.
* Update reference manual.
* Bugs fixed.

# ncProR 1.0.0

* Original release
