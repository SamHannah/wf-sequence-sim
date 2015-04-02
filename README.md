# wf-sequence-sim
Shiny app that runs simulations of a word-frequency sequence-effect experiment, based on my actual experiment; in the experiment and simulation, people/computer identifies a structured stream of high- and low-frequency words (HF & LF words) as either "old" (studied words, or "targets") or "new" (unstudied words, or "lures"). 

The simulation produces four graphs. First is a graph of the sequence effect  in terms of ∆P("old") (explained below); the second graph is a graph of the mean P("old") across LF lures, HF lures, HF targets and LF targets.  A mirror effect is indicated by bars that increase from left (LF lures) to right (LF targets). Error bars for both graphs are 95% confidence intervals. The third graph is a plotting of the criterion position (in Minerva-AL's intensity units; intensity being a proxy for familiarity) across trials for 2 single simulation, one in which the starting criterion position is low (.25 or lower, with mean intensity usually ~0.55 or slightly lower) or is high (.75 or higher).  The fourth graph shows the distribution of intensity units for each of the four test item types (LF lures, HF lures, HF targets, LF targets)

Simulation allows a person to set L, A, SF & dp parameters, and run 5, 25, 100 or 200 iterations/simulation. On reasonably fast machine, 200 iterations take just under an hour, five take about two to three minutes.  The simulation is set to default values of L = .75, A = .5, SF = 10, dp = 0.01, with 5 simulations.  This produces both a word-frequency mirror effect, and a reliable sequence effect in the first quartile only; this pattern is close to what is observed in the actual experiment. After setting desired parameter, iteration values, simulation is initiated by clicking on the "Run simulation" button on the bottom right of the sidebar.

Experiment details:
Study of 160-word list consisting of equal mix HF & LF words, then test with 160 targets and 160 lures, the latter also being an equal mix of HF & LF words.  Test words presented in quartets consisting of a pair of LF words followed by a pair of HF words.  For each pair, test status (target or lure) is constant.  So both LF pair members are targets, or both HF pair members are lures, etc.. HF pairs are crossed with LF pairs: LF(target)-HF(target); LF(target)-HF(lure); LF(lure)-HF(target); LFlure-HFtarget. There are 80 of each quartet type; quartet order is randomized at start of each session.  

A sequence effect is produced when HF words (either targets or lures) are more likely to be called "old" when preceeded by LF lures than when preceeded by LF targets, so the P("old") for the HF members of LF(lure)-HF(target/lure) quartets > P("old") for LF(lure)-HF(target/lure) quartets. Or, ∆P("old") > 0; ∆P("old") for HF(targets) = P("old)HF(target)|LF(lure) - ("old)HF(target)|LF(target).

The recognition test is broken into trial quartiles of 80 trials each, with five quartets of each type. Therefore, for HF test items, there are 10  stimuli of each Test Status (target/lure) X Prior Status (LF target/LF lure) condition occurring within each quartile.  For example, there are 10 HF lures that are preceeded by LF targets in each quartile, and there are also 10 HF lures that are preceeded by LF lures.

Simulation details:
The model is Minerva-AL (Jamieson, Hannah & Crump,2010; Jamieson, Crump & Hannah, 2012) supplemented by the criterion-calibration algorithm for shifting criterion on each trial (Hannah & Allan, 2011).  The Minerva-AL model has two parameters: L & A. Both control the proportion of information that is available from memory. L controls the proportion of information encoded at study. A controls the proportion of pre-study information stored in memory, and thus accessible at test.  Encoding uses Minerva-AL's discrepancy-encoding rule: Each probe presentation triggers a retrieval, which in turn produces an echo. The discrepancy between this echo and the study probe is what is encoded, probabilistically, based on L.

The criterion calibration algorithm has two parameters: SF and dp.  The SF parameter controls the overall lability of criterion shifting--how big a jump the criterion takes in response to a stimulus. The dp parameter changes the SF parameter by some amount on each trial. Typically, 0<= dp < 1, so dp reduces SF by some proportion on each trial. This implements the idea that with each trial people become a bit more confident about where criterion position should go, and thus less willing to shift.  When dp <= 1, SF asymptotes on 0; When dp > 1, SF alternates between a negative and positive values.  Whith dp < 2, SF still heads to 0; when dp = 2, SF alternates between +/- values of original value.

