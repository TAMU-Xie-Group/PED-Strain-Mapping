# PED Strain Mapping

### Abstract
In this work, we developed a method using the precession electron diffraction data to
map the residual elastic strain at the nano-scale. The diffraction pattern of each pixel was
first collected and denoised. Template matching was then applied using the center spot as
the mask to identify the positions of the diffraction disks. Statistics of distances between
the selected diffracted disks enable the user to determine “strain-free” baseline and to
generate strain maps. Strain mapping on an unstrained single crystal sapphire shows the
standard deviation of strain measurement is 0.5%. This approach does not require the user to select a
strain-free area as a reference and can work on datasets even the crystals oriented away
from zone axes. This method is expected to provide a robust way to study the residual
strain of various material systems that complements the existing algorithms for strain
mapping.


### Architecture
The strain mapping algorithm was written in Python and consists of five steps. 
First, it reads the .blo file using HyperSpy. Second, the
user can select a filter (Gaussian, or non-local means, or Wiener) to denoise the
diffraction patterns for each pixel. Third, the center diffraction disk is used as a template
to identify the positions of other diffracted beam disks using a correlation coefficient
cutoff (0.87 in the default setting). Fourth, the user will select a diffraction disk of
interest. The algorithm will go through all diffraction patterns, calculate the distance of
the selected disks, and generate a distance histogram. Note the distance is measured in
reciprocal space, and the unit is pixel. Finally, the user will determine a “strain-free
distance” based on the distance histogram, and the algorithm will create the “strain heat
map”.


### Instructions 
Download both Python files and run (each program has its own UI).
Video tutorial is available on YouTube: https://www.youtube.com/watch?v=0osP1SG77tY
