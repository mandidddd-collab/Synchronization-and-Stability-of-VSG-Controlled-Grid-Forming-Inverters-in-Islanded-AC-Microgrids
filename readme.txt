CASE-STUDY CODE
===============

This folder contains the MATLAB and Simulink code used for the case studies in the paper on synchronization and local stability of VSG-controlled islanded AC microgrids.

The code is divided into three folders.


1. lamuda
----------

This folder computes the critical loading factors for six benchmark systems:

  - WSCC 9-bus system
  - IEEE 14-bus system
  - IEEE 30-bus system
  - New England 39-bus system
  - IEEE 57-bus system
  - IEEE 118-bus system

It compares the existing and proposed synchronization conditions and exports
the numerical results and LaTeX tables used in the paper.

See lamuda/readme.txt for detailed instructions.


2. root_locus
-------------

This folder calculates and plots the dominant-pole root locus of the IEEE 39-bus system. The script follows the synchronous equilibrium as the loading factor changes, calculates the reduced-VSG eigenvalues, and identifies the critical loading point associated with the phase-angle boundary.

The folder includes the MATLAB script, numerical data, and generated
root-locus figures.

See root_locus/readme.txt for detailed instructions.


3. simulink
-----------

This folder contains the IEEE 39-bus Simulink model, with the source and load buses modeled as VSGs. It is used to reproduce the sequential-loading frequency-response case study in the paper.

The folder includes the Simulink model, the MATLAB plotting script, saved simulation data, and generated frequency-response figures.

See simulink/readme.txt for detailed instructions.


Software requirements
---------------------

MATLAB is required for the lamuda and root_locus folders.

MATLAB, Simulink, and Simscape Electrical Specialized Power Systems are required for the simulink folder.
