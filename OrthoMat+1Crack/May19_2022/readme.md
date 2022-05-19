## Description:

- The main file is _'meka_orthoSCA3D.for'_
- The formulation is present in _'meka_formulation.for'_
- files _'meka_standard_support.for'_ and _'meka_matrix_inverse.for'_ are additional files that contain standard subroutines that we need (like calculation of principal values and directions, calculation of inverse of a matrix, finding solutions for a system of linear equations etc.).

## Comments:
1. For Transverse Tension(Mode-3) and Compression(Mode-4) we assume the crack plane like we did for Longitudinal case, instead of finding the crack plane through the principal stresses and strains calculation
2. Shear23 case is not formulated in the code.  

## MODES:
- MODE = 1:   [Longitudinal Tension]
- MODE = 2:   [Longitudinal Compression]
- MODE = 3:   [Transverse Tension]
- MODE = 4:   [Transverse Compression]
- MODE = 7:   [Shear12]
- MODE = 8:   [Shear13]

- [x] **Failure Criteria**: Maximum Stress Failure criteria

## Shear13: Results

<img src="https://github.com/MekaSaiKrishna/CMC/blob/UMAT/OrthoMat%2B1Crack/May19_2022/Images/StressesVsTime_S13.svg">

<p float="left">
  <img src="https://github.com/MekaSaiKrishna/CMC/blob/UMAT/OrthoMat%2B1Crack/May19_2022/Images/S13vsSeparation_S13.svg" width="400" />
  <img src="https://github.com/MekaSaiKrishna/CMC/blob/UMAT/OrthoMat%2B1Crack/May19_2022/Images/S13vsepscr13_S13.svg" width="400" />
  <img src="https://github.com/MekaSaiKrishna/CMC/blob/UMAT/OrthoMat%2B1Crack/May19_2022/Images/S33vstime_S13.svg" width="400" /> 
  <img src="https://github.com/MekaSaiKrishna/CMC/blob/UMAT/OrthoMat%2B1Crack/May19_2022/Images/S11vstime_S13.svg" width="400" />
</p>

### Input Properties:
Material Property  | Values
------------- | -------------
E1      | 128 GPa
E2      | 7.6 GPa
G12     | 4.4 GPa
nu12    | 0.35
G23     | 2.62 GPa
GIC_F   | 40 KJ/m^2
GIIC_F  | 40 KJ/m^2
GIC_M   | 40 KJ/m^2
GIIC_M  | 40 KJ/m^2
XT      | 2.3 GPa
XC      | 1.531 GPa
YT      | 44 MPa
YC      | 44 MPa
TAUcr12 | 78.4 MPa
TAUcr23 | 78 MPa
damp    | 0
    


