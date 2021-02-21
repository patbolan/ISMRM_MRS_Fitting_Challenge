20210220
Patrick Bolan
This is my code for the ISMRM MRS Fitting Challenge of 2016. This code is relatively simple: each spectrum is modeled as a linear combination of the provided basis sets. All metabolites are modified by a zero-order phase, a delta-frequency, a Lorentzian broadening term, and a Gaussian broadening term (all global, not per-metabolite). Essentially this gives a Voigt lineshape rather than the more flexible model described by Provencher in LC Model. Note that  water lineshape correction was not implemented.

This code runs on Matlab2017b, and probably other versions as well (I did not record the version used for the official challenge submission). The optimizer is Matlab's lsqnonlin() function, and the goal was minimization of the real & imaginary residuals over a specified frequency range. For convenience I used a MRS toolset called MBS (Minnesota Breast Spectroscopy), written by me with guidance from Mike Garwood for all of our papers on breast MRS (starting with DOI: 10.1002/mrm.10654). This object-based library is included in the code. 

I ran the code one sample at a time using go_process_one_case.m, manually adjusting a few parameters (HSVD water suppression, include lipid, pre-phase the water peak) as needed. I did not record which settings I used for each case.

Regarding motivation: I submitted data for this challenge to see how well this very simple approach would perform compared to the more sophisticated and well-known tools. I did not implement baseline correction or lineshape correction (other than Voigt) to maintain that simplicity. But I should have! Might have been more competitive.

This code is public domain, but if you find it useful please cite our MRM 2003 paper, or include Mike and I in your acknowledgments. Please let me know if you have questions or comments: bola0035 at umn dot edu.

