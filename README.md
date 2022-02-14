# CMC
UMAT Codes related to my research

## Problem-1: Orthotropic - Elastic
Description: There are no cracks, it models the orthotropic material behaviour. 
* *matrixInverse.for* - responsible for inversion of matrices
* *myStandardSupport.for* - to find principal stresses and strains
* *ortho_elastic.for* - Subfile1 where all the governing equations are incorporated
* *Ortho3D* - the Main file which is used in the ABAQUS while submitting job

- state variables: 12
- Material Properties input order: E<sub>1</sub>, E<sub>2</sub>, E<sub>3</sub>, G<sub>12</sub>, G<sub>13</sub>, G<sub>23</sub>, G<sub>12</sub>, nu<sub>12</sub>, nu<sub>13</sub>, nu<sub>23</sub>

## Problem-2: Orthotropic SCA (Smeared Crack Approach)
Description: There is only "1" crack and there is no damping matrix i.e. D<sup>da</sup>=[0]<sub>3X3</sub> considered.
The crack plane's normal is along the '1' direction.
* *matrixInverse.for* - responsible for inversion of matrices
* *misesAndSmearedCrack_ortho.for* - Subfile1, where all the governing equations are incorporated
* *myStandardSupport_orthoSCA.for* - to find principal stresses and strains
* *orthoSCA3D.for* - the Main file which is used in ABAQUS while submitting job

- state variables:
- Material Properties input order: E<sub>1</sub>, E<sub>2</sub>, E<sub>3</sub>, G<sub>12</sub>, G<sub>13</sub>, G<sub>23</sub>, G<sub>12</sub>, nu<sub>12</sub>, nu<sub>13</sub>, nu<sub>23</sub>, sig<sup>cr0</sup>, tau<sub>1</sub><sup>cr0</sup>,tau<sub>2</sub><sup>cr0</sup>
 G<sub>IC</sub>, G<sub>IIC</sub>
