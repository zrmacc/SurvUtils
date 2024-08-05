## Version 0.9.0

	* Added functions for estimation and inference on cumulative incidence curves (CICs).
		- Added `GenCRData` to simulate competing risks data.
		- Added `OneSampleCIC` to estimate the cumulative incidence at a given time point.
	* Renamed `GetCurves`to `SurvCurves` and added an analogous function `CICurves` for CICs. 
	* Updated `GenPseudo` to generate psuedo-values for the cumulative incidence.
	* Updated `GenData` and `GenCRData` to return simplified output by default.

## Version 0.8.4

	* Added `RMSTInfluence` to calculate influence functions for RMST.
	* Added `GenPseudo` to generate pseudo-values for the RMST and survival probability. 
		- Note: influence function and pseudo-value calculation becomes slow for sample sizes exceeding 1e4.
