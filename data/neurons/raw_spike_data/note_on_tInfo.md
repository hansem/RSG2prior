# Quick question: the time between Set and Go in tInfo (i.e., code value %6 - %5) does not match the value of ’t’ (production interval). How do you define t? 

 Here's more info about tInfo. First of all, as long as the discrepancy is in the order of 20ms, it doesn't matter as all analyses are using 20ms bin and 40ms smoothing kernel (it matters if you measure latency of visual responses, for example).

 One thing I noticed from MWork code is that t (behavioral tp) is defined from set onset (to be exact, when WMorks send the command of refreshing monitor) to when the target is acquired (eye position enters the target window). However, %6 in tInfo is when eye position becomes out of the fixation window (i.e. when MWorks send an analog signal to Blackrock through NIDAQ). %5 is when the photodiode detects the flash presented together with the Set and send a signal to blackrock. 

 So the discrepancy can come from difference in 1) the begin time (set onset vs 5%) as %5 could be delayed by the refresh cycle+potential delay from photodiode to blackrock, 2) the end time as saccade duration is included in t.

 I only used the photodiode when any visual stimuli changed (e.g. fixation presented, target on, ready, set, and reward) to get a better timing than the NIDAQ analog (so photodiode was not used for %6). I believe timing estimate from photodiode is more accurate for these visual stimuli changes than NIDAQ (as we don't have to worry about refresh cycle). But regarding 'Go', this could differ. When I plot PSTHs time-locked to 'Go', I used %6. But an alternative is use %5+ t (behavioral tp). Anyway, I don't think this matters due to the smoothing. 

 Hope now it's clear. It'd be easier to explain with connection diagram of MWorks, monitor, NIDAQ, and blackrock.