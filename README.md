# MPS-Nonlinear-Optics
There are three main files of this project, `MPS_MPO_fix.py`, `TimeEvolutionMultiThread.py`, `MultipleSuperMode.py`. 
The `MPS_MPO_fix.py` contains the basic Matrix product state functions, including transfer a general state to MPS and general operator to MPO. It also provide one mode and two mode operator apply on MPS
There are three two mode operator functions: 
`two_mode()` function read in the MPS, operator and location and return the new MPS after the operator. 
`two_mode_o()` function read in the MPS, operator and location and modify the MPS to new MPS after the operator.
`two_mode_s()` function read in the only part of MPS near the location of operator, operator and location and return the part of new MPS after the operator.
There is a `MPSSample.ipynb` for showing how `MPS_MPO_fix.py` works with examples.


The `AttemptTimeEvolution.py` is the original python file for calculating the time evolution of MPS begin with displaced state of the supermode. And it will plot the photon density plot and save the MPS per 100 iteration.
The `TimeEvolutionMultiThread.py` is the multiple thread version of `AttemptTimeEvolution.py`, which could be faster or maybe not. With my own computer (cpu i7 10750 12 cores), the speed is about 5 to 10 times faster. However on the server panther2 the multiple thread version is slower.
And the time evolution of MPS will be save as `.p` (pickle file)


After we got the pickle file of MPS, we use the `MultipleSuperMode.py` to plot the Wigner function, calculate Negativity, purity, and energy.
And the `AttemptTimeEvolution.py` is the file that only plot one Wigner functions instead all, which is clear to read.
There are also `CheckDisplacementOperator.py` file for checking whether our wigner function is accurate for displaced state.


The `THDE.py` file is to plot the THDE approximation of the state.
