---
title: "Basis set convergence of molecular properties: Geometry and Vibrational Frequency"
teaching: 30
exercises: 60
questions:
- "Do the calculated molecular properties of a molecule converge with increasing basis set size?"
objectives:
- "Perform geometry optimization calculations."
- "Perform frequency calculations."
- "Plot the value of a molecular property vs. basis set."
keypoints:
- "It is only correct to perform a frequency calculation if you are at a local energy minimum conformation for your molecule."
- "The values of molecular properties *should* converge as you increase the size of the basis set."
---
## Overview
In this exercise we will perform geometry optimization and vibrational frequency analysis of a small molecule using three different basis sets of increasing size. Theoretically the molecular geometry and vibrational frequencies should converge as the size of the basis set increases. We’ll see if this is true!

For this exercise we will start by looking at the geometry and vibrational frequencies of a simple linear molecule, carbon dioxide (CO2).  This exercise is especially relevant because three of CO2’s four vibrational modes – or two of its three vibrational frequencies (because two of the modes have the same frequency) – are excitable by the infrared radiation emitted by Earth’s surface and this is what makes CO2 a potent greenhouse gas on our planet.  

In order for a vibrational mode to be excited by electromagnetic radiation, it must involve a change in the molecule’s dipole moment.  A visualization of CO2’s vibrational modes is available here: http://www.chemtube3d.com/vibrationsco2/.  Just by looking at the animations of the vibrational modes: which modes cause a change in the dipole moment of CO2?  Do any of them have identical frequencies?

> ## Exercise
> Construct a CO<sub>2</sub> molecule.  You can use a molecular visualization program like Avogadro, create your own set of x,y,z coordinates, or use a z-matrix specification.
>
>> ## Solution
>> One possible set of x,y,z coordinates
>> ~~~
>> O       -1.25000        0.00000        0.00000
>> C        0.00000        0.00000        0.00000
>> O        1.25000        0.00000        0.00000
>> ~~~
>> {: .output}
>>
>> One possible z-matrix
>> ~~~
>> C
>> X 1 1.0
>> O 1 1.25 2 90.0
>> 0 1 1.25 2 90.0 3 180.0
>> ~~~
>> {: .output}
> {: .solution}
{: .challenge}

## Calculation setup

First we import the required python functions into the notebook.
~~~
import psi4
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
~~~
{: .language-python}

Next, set up the molecule for your calculation. Don't forget to include your molecular charge and multiplicity!
~~~
psi4.set_memory('2 GB')
psi4.set_num_threads(2)

co2 = psi4.geometry("""
symmetry c1
0 1
YOUR MOLECULE SPECIFICATION GOES HERE
""")
~~~
{: .language-python}

## Geometry optimization and vibrational frequency analysis
Optimize the geometry of the molecule to ensure it is at its minimum energy conformation.
~~~
psi4.set_output_file('co2_geomopt_HF-DZ.dat', False)
psi4.set_options({'g_convergence': 'gau_tight'}) # this forces the optimizer to get close to the minimum
psi4.optimize('scf/cc-pVDZ', molecule=co2, dertype='gradient')
~~~
{: .language-python}
~~~
Optimizer: Optimization complete!
-187.65250937052735
~~~
{: .output}

Let's look at the molecular geometry after optimization and compute the C-O bond length "by hand", i.e., taking {x,y,z} coordinates for one of the O atoms and the C atom and computing the Euclidean distance between them.
~~~
scf_dz_xyz = 0.529177*np.array(co2.geometry()) # the factor of 0.529177 converts from Bohr to Ångstrom
print('The atomic coordinates after optimization are:\n', scf_dz_xyz)

scf_dz_bondlen = np.linalg.norm(scf_dz_xyz[1,:]-scf_dz_xyz[0,:]) # look up NumPy's linalg.norm() function!
print('The C-O bond length after optimization is:\n', scf_dz_bondlen)
~~~
{: .language-python}

~~~
The atomic coordinates after optimization are:
 [[-1.14055067e+00  6.48113657e-17  1.62668842e-15]
 [-1.13961561e-11 -1.72775378e-16 -4.33645706e-15]
 [ 1.14055067e+00  6.48113666e-17  1.62668842e-15]]
The C-O bond length after optimization is:
 1.1405506737107873
 ~~~
 {: .output}

The command below will execute a vibrational frequency analysis (requiring calculation of the Hessian) for this molecule. We will request that Psi4 approximate the Hessian by calculating finite differences of the gradient. This method is often faster than directly calculating the Hessian, and second derivatives of the energy are often unavailable for various quantum chemistry methods.
~~~
psi4.set_output_file('co2_vibfreq_HF-DZ.dat', False)
scf_dz_energy, scf_dz_wfn = psi4.frequency('scf/cc-pVDZ', molecule=co2, return_wfn=True, dertype='gradient')
~~~
{: .language-python}

~~~
9 displacements needed.
 1 2 3 4 5 6 7 8 9
~~~
{: .output}

For a linear molecule, we expect there to be 3N-5 normal modes (and therefore normal mode frequencies). In the case of CO 2 , we expect (3x3) - 5 = 4 normal modes. We will print the computed vibrational frequencies out below.
~~~
for i in range(4):
    print(scf_dz_wfn.frequencies().get(0,i))
~~~
{: .language-python}

~~~
761.2953702942335
761.2953703599782
1513.3176259835404
2580.2533937232515
~~~
{: .output}

It looks like the first two frequencies are identical (as suggested by our resource: http://www.chemtube3d.com/vibrationsCO2.htm). Let's exclude the first one from the list of frequencies that we'll store.
~~~
scf_dz_vibfreq = [] # make an empty list to store frequencies

for i in range(3):
    scf_dz_vibfreq.append(scf_dz_wfn.frequencies().get(0,i+1))

print(scf_dz_vibfreq)
~~~
{: .language-python}

~~~
[761.2953703599782, 1513.3176259835404, 2580.2533937232515]
~~~
{: .output}

> ## Exercise: Vibrational Frequency Analysis with Larger Basis sets
> Using the commands you learned above, calculate the vibrational frequencies of CO<sub>2</sub> using the cc-pVTZ and cc-pVQZ basis sets.
>
>> ## Solution
>> For the cc-pVTZ basis set
>> ~~~
>> psi4.set_output_file('co2_geomopt_HF-TZ.dat', False)
>> psi4.optimize('scf/cc-pVTZ', molecule=co2, dertype='gradient')
>>
>> scf_tz_xyz = 0.529177*np.array(co2.geometry())
>> scf_tz_bondlen = np.linalg.norm(scf_tz_xyz[1,:]-scf_tz_xyz[0,:])
>>
>> psi4.set_output_file('co2_vibfreq_HF-TZ.dat', False)
>> scf_tz_energy, scf_tz_wfn = psi4.frequency('scf/cc-pVTZ', molecule=co2, return_wfn=True, dertype='gradient')
>>
>> scf_tz_vibfreq = []
>> for i in range(3):
>>     scf_tz_vibfreq.append(scf_tz_wfn.frequencies().get(0,i+1))
>> ~~~
>> {: .language-python}
>>
>> ~~~
>> Optimizer: Optimization complete!
>>  9 displacements needed.
>>  1 2 3 4 5 6 7 8 9
>> ~~~
>> {: .output}
>>
>> For the cc-pVQZ basis set
>> ~~~
>> # for cc-pVQZ
>> psi4.set_output_file('co2_geomopt_HF-QZ.dat', False)
>> psi4.optimize('scf/cc-pVQZ', molecule=co2, dertype='gradient')
>>
>> scf_qz_xyz = 0.529177*np.array(co2.geometry())
>> scf_qz_bondlen = np.linalg.norm(scf_qz_xyz[1,:]-scf_qz_xyz[0,:])
>>
>> psi4.set_output_file('co2_vibfreq_HF-QZ.dat', False)
>> scf_qz_energy, scf_qz_wfn = psi4.frequency('scf/cc-pVQZ', molecule=co2, return_wfn=True, dertype='gradient')
>>
>> scf_qz_vibfreq = []
>> for i in range(3):
>>     scf_qz_vibfreq.append(scf_qz_wfn.frequencies().get(0,i+1))
>> ~~~
>> {: .language-python}
>>
>> ~~~
>> Optimizer: Optimization complete!
>>  9 displacements needed.
>>  1 2 3 4 5 6 7 8 9
>> ~~~
>> {: .output}
> {: .solution}
{: .challenge}

> ## Exercise: Plotting the results
> Plot the C-O bond length and the vibrational frequencies vs. basis set size.
>
>> ## Solution
>> The C-O bond length
>> ~~~
>> plt.plot([1, 2, 3], [scf_dz_bondlen, scf_tz_bondlen, scf_qz_bondlen],
>>          marker='o', color='RebeccaPurple')
>> plt.xticks([1, 2, 3], ['DZ', 'TZ', 'QZ'])
>> plt.xlabel('Basis set')
>> plt.ylabel(r'Bond distance ($\AA$)')
>> plt.title('C-O bond distance in carbon dioxide')
>> plt.show()
>> ~~~
>> {: .language-python}
>> <img src="../fig/CO-bond.png"/>
>>
>> The frequencies
>> ~~~
>> # plot vibrational frequencies
>> plt.plot([1, 2, 3], [scf_dz_vibfreq[0], scf_tz_vibfreq[0], >> scf_qz_vibfreq[0]],
>>          marker='o', color='red')
>> plt.xticks([1, 2, 3], ['DZ', 'TZ', 'QZ'])
>> plt.xlabel('Basis set')
>> plt.ylabel(r'Frequency (cm$^{-1}$)')
>> plt.title(r'CO2 HF vibrational frequency #1')
>> plt.show()
>> ~~~
>> {: .language-python}
>> <img src="../fig/freq1.png"/>
>>
>> ~~~
>> plt.plot([1, 2, 3], [scf_dz_vibfreq[1], scf_tz_vibfreq[1], scf_qz_vibfreq[1]],
>>          marker='o', color='green')
>> plt.xticks([1, 2, 3], ['DZ', 'TZ', 'QZ'])
>> plt.xlabel('Basis set')
>> plt.ylabel(r'Frequency (cm$^{-1}$)')
>> plt.title(r'CO2 HF vibrational frequency #2')
>> plt.show()
>> ~~~
>> {: .language-python}
>> <img src="../fig/freq2.png"/>
>>
>> ~~~
>> plt.plot([1, 2, 3], [scf_dz_vibfreq[2], scf_tz_vibfreq[2], scf_qz_vibfreq[2]],
>>          marker='o', color='blue')
>> plt.xticks([1, 2, 3], ['DZ', 'TZ', 'QZ'])
>> plt.xlabel('Basis set')
>> plt.ylabel(r'Frequency (cm$^{-1}$)')
>> plt.title(r'CO2 HF vibrational frequency #3')
>> plt.show()
>> ~~~
>> <img src="../fig/freq3.png"/>
>>
> {: .solution}
{: .challenge}

> ## Comparing your results to experiment
> To compare your results to experiment, you could create a plot of the percent errors in the vibrational frequencies with respect to experiment.  Note that this will require you to use the experimental data and scale factors from the NIST CCCBDB web pages listed above.  Here is some starting code.
>
> ~~~
> scaled_vibfreq1 = np.array([0.908*scf_dz_vibfreq[0], 0.910*scf_tz_vibfreq[0], 0.908*scf_qz_vibfreq[0]])
> ### two more lines with new variables for the other two frequencies ###
>
> plt.plot([1, 2, 3], 100*(scaled_vibfreq1 - 667)/667, label='Vib. Freq. #1', color='red')
> ### two more lines with calls to plt.plot for the other two frequencies ###
> plt.xticks([1, 2, 3], ['DZ', 'TZ', 'QZ'])
> plt.xlabel('Basis set')
> plt.ylabel('Percent error')
> plt.title('Percent error in scaled HF vibrational frequencies')
> plt.legend()
> plt.show()
> ~~~
>{: language-python}
{: .callout}

## Project Extension: Geometry Optimization and Frequency analysis with post-Hartree Fock Methods
You can now repeat this procedure to with DFT or MP2 to study the effects of electron correlation on your convergence.  To do this, change all instances of scf to b3ylp (a very commonly used DFT functional) or mp2 (perhaps the most commonly used post-HF method).  Either one of these is a reasonable first approach to accounting for the effects of electron correlation.

Does DFT/B3LYP or MP2 generate frequencies that are in better agreement with experiment than HF?  To answer this question you’ll need to find the appropriate multiplicative scale factors to use with these methods.
