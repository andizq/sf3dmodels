Modelling multiple star forming regions
=======================================

**Example 1.** I will use the last two examples to illustrate how to join them in a *global grid*. The spatial region that is shared by two or more *sub-models* will inherit physical properties by weighting them with the local density, as explained in section 3.2 of Izquierdo et al (2018).

**The execution codes for both star forming regions are identical until just before the "writing" section.**

As there is no longer a single model, geometric changes will probably be required in each sub-model to better reproduce real scenarios. Let's add a couple of lines in the latter codes to account for the centering, inclination and systemic velocity of each region.

In the first example:
