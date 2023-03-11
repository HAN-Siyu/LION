# 0.2.9
Fix bug in `evaluatePrediction()`

# LION 0.2.8
mtry can be tuned.

# LION 0.2.7
Update models.

# LION 0.2.6

* Update reference manual.
* Bugs fixed.

# LION 0.2.5

* Update reference manual.
* Bugs fixed.

# LION 0.2.4

* Users can use `run_confidentPrediction()` to get prediction results of multiple methods in one run.
* Update reference manual.
* Bugs fixed.

# LION 0.2.3

* Add `run_LncADeep()` for prediction using LncADeep's feature set.
* Add `computeMLC()` to compute the most-like CDS (MLC) region.
* Add parameter `EDP` which supports applying entropy density profile to sequence intrinsic features.
* Update reference manual.
* Bugs fixed.

# LION 0.2.2

* Matthews correlation coefficient (MCC) is included into evaluation metrics.
* TP, TN, FP and FN are reported in evaluation results.
* Add `evaluatePrediction()` for prediction evaluation.
* Bugs fixed.

# LION 0.2.1

* Add argument verification.
* Update reference manual.
* Bugs fixed.

# LION 0.2.0

* `"retrain"` and `"feature"` modes are added to function `run_LION()`, `run_RPISeq()`, `run_lncPro()` and `run_rpiCOOL()`. Now users can use ready-to-use functions to extract features or retrain the models of LION, RPISeq, lncPro and rpiCOOL.
* External parallel cores can be passed to functions.
* Internal changes to improve efficiency.
* Update reference manual.
* Bugs fixed.

# LION 0.1.0

* Original release
