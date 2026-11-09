scripts in this folder and their functions:
- iptw: calculates weights for observations given a treatment model, takes as input a baseline dataset
- boot_iptw: does the same as iptw, takes as input longitudinal dataset and then slices to baseline to compute weights, given a treatment model. Gives as output only the ate of the outcome model specified, can take as input tte and mixed models both
- gcomp: performs g-computation given an outcome model, can take as input both longitudinal or baseline dataset, and performs predictions for tte and mixed models both
- interpolation: given a dataset and a time in point at which we wish to cut, performs interpolation, to use as many data points as possible
- alignment_psm_matching: performs matching in a similar fashion as matchit from the MatchIt library, but considering also time alignment, matching techniques allowed are: greedy, optimal and ilp_optimal
- alignment_psm_matching_interpolation: does the same as alignment_psm_matching, but for longitudinal data
- exact_matching: performs exact matching in a similar fashion as matchit from the MatchIt library, but considering also time alignment, matching techniques allowed are: greedy, optimal and ilp_optimal
- exact_matching_interpolation: performs exact matching, similarly to what has been done for the GBZ in VWM SAT 
- ilp_optimal: optimizer to perform optimal search given a caliper (option not allowed in matchit), while not throwing an error if the optimization step fails 
- boot_matching: wraps all possible matching for both endpoints and gives back the estimate of interest

In matching with duration for continuous events, the output can be complicated:
	- original: original dataset
	- aligned: takes the matched pairs and aligns them in the disease duration timeline
	- filtered: takes `aligned` as basis and filters time variable to keep only data after time zero, so for control the first observation may be months/years later
	- interpolated: takes `filtered` and adds the value of `y`, where it would have been, had it been measured for historical controls at baseline


To align correctly:
	- original: hui_linear_trt vs timen
	- aligned:  hui_linear_trt vs time_aligned (group by trt to see change in disease course)
	- filtered: hui_linear_trt vs time_aligned (can group by either trt or trt_overall)
	- interpolated: y vs time_aligned (can group by either trt or trt_overall)