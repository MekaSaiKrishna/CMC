# Documentation: 

The UMATs are for the case of an Orthotropic Material with 3 isotropic cracks (i.e. the D<sup>cr</sup> matrix is of Isotropic form)

* *matrixInverse.for* : computes inverse of a matrix, solves linear system of equations (*solveLinSysLU()*)
* *myStandardSupport.for* : computes principal stress values/principal directions
* *misesAndSmearedCrack_ortho.for* : SCA for Ortho (transversely isotropic to be precise) material with 3 iso cracks. D<sup>da</sup> ≠ 0, it is computed by assuming a damping term. The 3 crack planes are mutually perpendicular with their normals along the global co-ordinate axes.

  * Inputs in 'Material Properties' in CAE file:
    * E1, E2, G12, nu<sub>12</sub>, nu<sub>23</sub>, sigma0<sub>cr</sub>, GIC, damp
