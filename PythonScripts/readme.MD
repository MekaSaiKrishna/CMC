# README file about the Python Scripts of ABAQUS


## *tiedLaminate1.py*: 
**Description**: 4 ply laminate with tied interaction between plies and aligned meshing [partitioning]. 
The partitioning is done by giving few parameters as the input, the layup sequence is [+45/-45/-45/+45].

**Limitations**
- When thickness of the strip is 1mm or lower the assignment of sections is fucked
- Material Properties: Both Lamina and Matrix are SAME!!
- interactions: [TIE CONSTRANT]
- The code for the partitionining can be improved

**Further Improvements to do:**
- Constraints: Cohesive Zone to be introduced
- Material Properties: Lamina and Matrix to be differentiated
- Meshing: Transition region and assign different meshing algos to different sections

## *CZMlaminate1.py*:
**Description**: 

**Limitations**
- j
- k

**Further Improvements to do:**

The UMAT file for SCA used in these scripts are here. [^1]

[^1]: These are the SCA files that ABAQUS read: [link](https://github.com/MekaSaiKrishna/CMC/tree/UMAT/MekaSCAcodes).
