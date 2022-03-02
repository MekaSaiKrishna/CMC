# Documentation: 

The UMATs are for the case of an Orthotropic Material with 3 Orthotropic cracks (i.e. the D<sup>cr</sup> matrix is of Orthotropic form)

* *matrixInverse.for* : computes inverse of a matrix, solves linear system of equations (*solveLinSysLU()*)
* *myStandardSupport.for* : computes principal stress values/principal directions
* *misesAndSmearedCrack_ortho.for* : SCA for Ortho (transversely isotropic to be precise) material with 3 ortho cracks. D<sup>da</sup> ≠ 0, it is computed by assuming a damping term. The 3 crack planes are mutually perpendicular with their normals along the global co-ordinate axes.
* *USCA3D.for* : main subroutine file which is fed to the CAE file.

  * Inputs in 'Material Properties' in CAE file:
    * E1, E2, G12, nu<sub>12</sub>, nu<sub>23</sub>, sigma0<sub>cr</sub>, GIC, GIIC, tau1<sub>cr</sub>0, tau2<sub>cr</sub>0, damp
  * Element size: *1mm*     



Displacement boundary conditions are used on the faces [^1].

[^1]: U1=0 for a face with normal along '1' and similarly U2 and U3 are 0 for one of the two faces having normals along direction 2 and 3. Another face which has a normal along '1' is given a displacement U1=0.06. So this is uniaxial tensile loading along '1' direction on this element.
