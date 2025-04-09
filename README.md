# EFMouse: a Matlab toolbox to model stimulation-induced electric fields in the mouse brain
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14777232.svg)](https://doi.org/10.5281/zenodo.14777232) ![License](https://img.shields.io/badge/License-MIT-yellow.svg)

The [EFMouse.m](https://github.com/klabhub/EFMouse/blob/main/EFMouse.m) Matlab class models electric fields induced by current stimulation in
the mouse.

For details see [Sanchez-Romero et al. (2024). bioRxiv.](https://doi.org/10.1101/2024.07.25.605227) <br /> 
Run the Matlab notebooks [montage_4x1.mlx](https://github.com/klabhub/EFMouse/blob/main/montage_4x1.mlx) and [montage_1x1.mlx](https://github.com/klabhub/EFMouse/blob/main/montage_1x1.mlx) to reproduce results in Sanchez-Romero et al.<br />
For convenience, you can visualize (not interact with) the notebooks directly here [montage_4x1.hmtl](https://klabhub.github.io/EFMouse/montage_4x1.html) and [montage_1x1.hmtl](https://klabhub.github.io/EFMouse/montage_1x1.html).

Developed by Ruben Sanchez-Romero, Sibel Akyuz and Bart Krekelberg<br /> 
Center for Molecular and Behavioral Neuroscience (CMBN), Rutgers Newark<br/> 

For support open an [issue](https://github.com/klabhub/EFMouse/issues).

If you use EFMouse in your research, please cite our manuscript: Sanchez-Romero R., Akyuz, S., & Krekelberg, B. (2024). EFMouse: a Matlab toolbox to model stimulation-induced electric fields in the mouse brain. bioRxiv. https[]()://doi.org/10.1101/2024.07.25.605227

## Installation and dependencies
Clone the repository from github `git clone https://github.com/klabhub/EFMouse.git`

Analyzing results with reference to the Allen mouse atlas requires the open-source stand-alone [FSL 
package](https://fsl.fmrib.ox.ac.uk/fsl/docs/#/). We recommend using FSL 6.7.0 onwards. (FSL applies a transformation to the simulation results to map them to the Allen atlas.)

Mesh operations require Matlab Partial Differential Equation Toolbox.

## Tutorial
To get started, open [montage_4x1.mlx](https://github.com/klabhub/EFMouse/blob/main/montage_4x1.mlx) or [montage_1x1.mlx](https://github.com/klabhub/EFMouse/blob/main/montage_1x1.mlx) in Matlab and follow the compute pipeline, the analysis suggestions, 
and for additional options look up the help for each of the functions used in the tutorial.
