# README file about the SCA UMAT 


### *meka_formulation.for*: 
**Description**: The formulation is done in this file

**Features**:
- Single crack in an orthotropic material, crack orientation is dependent on the mode of failure
- Failure Criteria: Maximum Stress Criteria

**Limitations**
- include characteristic length dependent on angle, not on Abaqus [not sure what this means]
- compression failure
- crack closing in compression
- For Transverse Tension(Mode-3) and Compression(Mode-4) we assume the crack plane like we did for Longitudinal case
- Shear23 case is **"NOT"** present in this code

**Further Improvements to do:**
- Try to capture the S22,S33,S23 failure using the method used by Haotian
- In the long run, formulate SCA equations by long deformation concepts
- Find a small bug in the code which leads to an error for the case of laminates where all the stresses start depreciating once one of the stresses reaches its peak


* *meka_orthoSCA3D.for*: The main file that ABAQUS reads
* *meka_matrix_inverse.for* : computes inverse of a matrix, solves linear system of equations (*solveLinSysLU()*)
* *meka_standard_support.for* : computes principal stress values/principal directions


  * Inputs in 'Material Properties' in CAE file:
    * E1, E2, G12, nu<sub>12</sub>, nu<sub>23</sub>, sigma0<sub>cr</sub>, GIC, damp
  * Element size: *1mm*     

### Property Variables List [ABAQUS UNITS]
| Variable      | Value |
| ----------- | ----------- |
| E1       |         128e3 |
| E2       |         7.6e3 |
| G12      |         4.4e3 |
| nu12     |          0.35 |
| G23      |        2.62e3 |
| GIC_F    |            40 |
| GIIC_F   |             4 |
| GIC_M    |             2 |
| GIIC_M   |             1 |
| XT       |         2.3e3 |
| XC       |       1.531e3 |
| YT       |            44 |
| YC       |            44 |
| TAUcr12  |          78.4 |
| TAUcr23  |            78 |
| damp     |             0 |


Displacement boundary conditions are used on the faces [^1].

[^1]: U1=0 for a face with normal along '1' and similarly U2 and U3 are 0 for one of the two faces having normals along direction 2 and 3. Another face which has a normal along '1' is given a displacement U1=0.06. So this is uniaxial tensile loading along '1' direction on this element.

[^1]: These are the SCA files that ABAQUS read: [link](https://github.com/MekaSaiKrishna/CMC/tree/UMAT/MekaSCAcodes).
