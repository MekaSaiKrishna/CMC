# Attempt 3

# Date Modified: June 22nd, 2022

# * Limitations: When thickness of the strip is 1mm or lower the assignment of sections is fucked
# * Material Properties: Both Lamina and Matrix are SAME!!
# * interactions: [TIE CONSTRANT]
# * Partitioning Plies is done (but can be improved)

# * Further Improvements to be done
# --> Constraints: Cohesive Zone to be introduced
# --> Material Properties: Lamina and Matrix to be differentiated
# --> Meshing: Transition region and assign different meshing algos to different sections

# NOTE: Note the path of the UMAT file and the path of the file under which the CAE file is going to get saved.

from abaqus import *
from abaqusConstants import *
import __main__

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

import numpy as np
import math

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

#------------------------------------------------------------------------------#
############____________INPUT AREA__________####################################

Seq = [45,-45,-45,45] # Layup Sequence
t   = 5               # thickness of matrix-cracking strips
w   = 0.25            # thickness of each ply
L   = 200.0           # side length of square laminate
N   = 10              # partition points on each side (excluding the extreme corners)
meshSize=14.0

numLayers = len(Seq)  #number of layers in laminate 
################################################################################

# Creating a Lamina Model

mdb.models.changeKey(fromName='Model-1', toName='Lamina1')
laminaModel=mdb.models['Lamina1']

#------------------------------------------------------------------------------
# Creating the part

# a) Sketch the plate using the rectangle tool
plateProfileSketch = laminaModel.ConstrainedSketch(name='Lamina Sketch',sheetSize=200.0)
plateProfileSketch.rectangle(point1=(-L/2,-L/2), point2=(L/2,L/2))

# b) Create 'Ply-j' where 'j' is the layer number
for i in range(numLayers):
    laminaPart = laminaModel.Part(dimensionality=THREE_D, name='Ply-'+str(i+1), 
        type=DEFORMABLE_BODY)
    laminaPart.BaseSolidExtrude(depth=w, sketch=plateProfileSketch)

#------------------------------------------------------------------------------
# Creating Local and Global DatumCsys for each Ply

for i in range(numLayers):
    GlobalCoordSys = laminaModel.parts['Ply-'+str(i+1)].DatumCsysByThreePoints(coordSysType=
        CARTESIAN, name='Global', origin=(0.0,0.0,0.0), point1=(L/2,0.0,0.0),
        point2=(L/2,L/2,0))

layupSeq=np.array(Seq)*np.pi/180

for i in range(numLayers):
    if Seq[i] > 0:
        point1_temp = ( (L/2)*np.cos(layupSeq[i]),(L/2)*np.sin(layupSeq[i]),w)
        point2_temp = (-(L/2)*np.cos(layupSeq[i]),(L/2)*np.sin(layupSeq[i]),w)
        LocalCoordSys = laminaModel.parts['Ply-'+str(i+1)].DatumCsysByThreePoints(coordSysType=
            CARTESIAN, name='Local', origin=(0.0,0.0,0.0), point1=point1_temp, point2=point2_temp)
    elif Seq[i] < 0:
        point1_temp = ((L/2)*np.cos(layupSeq[i]),(L/2)*np.sin(layupSeq[i]),w)
        point2_temp = ((L/2)*np.cos(layupSeq[i]),-(L/2)*np.sin(layupSeq[i]),w)
        LocalCoordSys = laminaModel.parts['Ply-'+str(i+1)].DatumCsysByThreePoints(coordSysType=
            CARTESIAN, name='Local', origin=(0.0,0.0,0.0), point1=point1_temp, point2=point2_temp)
    else:
        LocalCoordSys = laminaModel.parts['Ply-'+str(i+1)].DatumCsysByThreePoints(coordSysType=
        CARTESIAN, name='Local', origin=(0.0,0.0,0.0), point1=(L/2,0.0,0.0),
        point2=(L/2,L/2,0))

#------------------------------------------------------------------------------
# Creating Material and assigning material properties

# a) Defining 'LAMINA' material
laminaMaterial = laminaModel.Material(name = 'LAMINA')
laminaMaterial.Density(table=((1.59e-09, ), ))
laminaMaterial.Depvar(n=40)
laminaMaterial.UserMaterial(mechanicalConstants=(
   128000.0, 7600.0, 4400.0, 0.35, 2620.0, 40.0, 4.0, 2.0, 1.0, 2300.0, 
   1531.0, 44.0, 44.0, 78.4, 78.0,0.0))
laminaMaterial.userMaterial.setValues(mechanicalConstants=(
   128000.0, 7600.0, 4400.0, 0.35, 2620.0, 40.0, 4.0, 2.0, 1.0, 2300.0, 
   1531.0, 44.0, 44.0, 78.4, 78.0,0.0))

# b) Defining 'MATRIX' material
matrixMaterial = laminaModel.Material(name = 'MATRIX')
matrixMaterial.Density(table=((1.59e-09, ), ))
matrixMaterial.Depvar(n=40)
matrixMaterial.UserMaterial(mechanicalConstants=(
   128000.0, 7600.0, 4400.0, 0.35, 2620.0, 40.0, 4.0, 2.0, 1.0, 2300.0, 
   1531.0, 44.0, 44.0, 78.4, 78.0,0.0))
matrixMaterial.userMaterial.setValues(mechanicalConstants=(
   128000.0, 7600.0, 4400.0, 0.35, 2620.0, 40.0, 4.0, 2.0, 1.0, 2300.0, 
   1531.0, 44.0, 44.0, 78.4, 78.0,0.0))

#------------------------------------------------------------------------------
# Creating Sections

Lamina_Section = laminaModel.HomogeneousSolidSection(material='LAMINA', name=
   'LaminaSection',thickness=None)
Matrix_Section = laminaModel.HomogeneousSolidSection(material='MATRIX', name=
   'MatrixSection',thickness=None)

# !NOTE: Assigning sections is done after partitioning!
#------------------------------------------------------------------------------
# Create the Assembly

laminaAssembly = laminaModel.rootAssembly
for i in range(numLayers):
    laminaInstance = laminaAssembly.Instance(name='Ply-'+str(i+1)+'-inst', 
        part=laminaModel.parts['Ply-'+str(i+1)], dependent=ON)
    laminaModel.rootAssembly.instances['Ply-'+str(i+1)+'-inst'].translate(
        vector=(0.0,0.0,w*i))

#------------------------------------------------------------------------------
# Defining sets for boundary conditions

# Case-(1): Loading in X-direction 

# Loading Faces [C D G F]
pt1 = (L/2.,0.0,1*(w/2.))
init = laminaModel.rootAssembly.instances['Ply-1-inst'].faces.findAt((pt1,))
for i in range(numLayers):
    z = ((2*i)+1)*(w/2.)
    pt = (L/2.,0.0,z)
    loadSurf = laminaModel.rootAssembly.instances['Ply-'+str(i+1)+'-inst'].faces.findAt((pt,))
    if i==0:
        tot = init
    else:
        tot = tot+ loadSurf

BoundSet_CDGF = laminaModel.rootAssembly.Set(faces=tot, name='BoundSet_CDGF')

# Constrained Face [B A H E]
pt1 = (-L/2.,0.0,1*(w/2.))
init = laminaModel.rootAssembly.instances['Ply-1-inst'].faces.findAt((pt1,))
for i in range(numLayers):
    z = ((2*i)+1)*(w/2.)
    pt = (-L/2.,0.0,z)
    loadSurf = laminaModel.rootAssembly.instances['Ply-'+str(i+1)+'-inst'].faces.findAt((pt,))
    if i==0:
        tot = init
    else:
        tot = tot+ loadSurf

BoundSet_BAHE = laminaModel.rootAssembly.Set(faces=tot, name='BoundSet_BAHE')

# Fixed point [A]
ptA = (-L/2.,L/2.,0)
SetPt_A = laminaModel.rootAssembly.Set(name='BoundSet_A', vertices=
    laminaModel.rootAssembly.instances['Ply-1-inst'].vertices.findAt((ptA,),))

# Fixed point [B]
ptB = (-L/2.,-L/2.,0)
SetPt_B = laminaModel.rootAssembly.Set(name='BoundSet_B', vertices=
    laminaModel.rootAssembly.instances['Ply-1-inst'].vertices.findAt((ptB,),))
#------------------------------------------------------------------------------
# Partitioning Plies based on the angle of fiber orientation

d=L/(N+1)

leftTopCorner=((-L/2)*0.99,(L/2)*0.99,w)
rightBottomCorner=((L/2)*0.99, (-L/2)*0.99,w)

rightTopCorner=((L/2)*0.99,(L/2)*0.99,w)
leftBottomCorner=((-L/2)*0.99, (-L/2)*0.99,w)

# !NOTE!: Need to find the reason behind using datums[3]

#------------------------------------------------------------------------------
# Partitioning Plies based on the angle of fiber orientation

d=L/(N+1)

leftTopCorner=((-L/2)*0.99,(L/2)*0.99,w)
rightBottomCorner=((L/2)*0.99, (-L/2)*0.99,w)

rightTopCorner=((L/2)*0.99,(L/2)*0.99,w)
leftBottomCorner=((-L/2)*0.99, (-L/2)*0.99,w)

# !NOTE!: Need to find the reason behind using datums[3]
for j in range(numLayers):
    laminaPart = laminaModel.parts['Ply-'+str(j+1)]
    if Seq[j] > 0:
        # laminaPart = laminaModel.parts['Ply-'+str(j+1)]
        for i in range(N+1):
            laminaCells = laminaPart.cells
            selectCell  = laminaCells.findAt(rightBottomCorner)
            laminaPart.PartitionCellByPlanePointNormal(cells=
                selectCell,normal=laminaPart.datums[3].axis2,point=(((-L/2)+(i*d)),(-L/2)*0.99,w))
            laminaCells = laminaPart.cells
            selectCell  = laminaCells.findAt(rightBottomCorner)
            laminaPart.PartitionCellByPlanePointNormal(cells=
                selectCell,normal=laminaPart.datums[3].axis2,point=(((-L/2)+(i*d)+t),(-L/2)*0.99,w))
        
        for i in range(N):
           laminaCells = laminaPart.cells
           selectCell  = laminaCells.findAt(leftTopCorner)
           laminaPart.PartitionCellByPlanePointNormal(cells=
               selectCell,normal=laminaPart.datums[3].axis2,point=((-L/2)*0.99,((-L/2)+((i+1)*d)),w))
           laminaCells = laminaPart.cells
           selectCell  = laminaCells.findAt(leftTopCorner)
           laminaPart.PartitionCellByPlanePointNormal(cells=
               selectCell,normal=laminaPart.datums[3].axis2,point=((-L/2)*0.99,((-L/2)+((i+1)*d)+t),w))
    elif Seq[j] < 0:
        for i in range(N+1):
            laminaCells = laminaPart.cells
            selectCell  = laminaCells.findAt(rightTopCorner)
            laminaPart.PartitionCellByPlanePointNormal(cells=
                selectCell,normal=laminaPart.datums[3].axis2,point=(((-L/2)+(i*d)),(L/2)*0.99,w))
            laminaCells = laminaPart.cells
            selectCell  = laminaCells.findAt(rightTopCorner)
            laminaPart.PartitionCellByPlanePointNormal(cells=
                selectCell,normal=laminaPart.datums[3].axis2,point=(((-L/2)+(i*d)+t),(L/2)*0.99,w))
        for i in range(N):
            laminaCells = laminaPart.cells
            selectCell  = laminaCells.findAt(leftBottomCorner)
            laminaPart.PartitionCellByPlanePointNormal(cells=
                selectCell,normal=laminaPart.datums[3].axis2,point=((-L/2)*0.99,((L/2)-(i+1)*d),w))
            laminaCells = laminaPart.cells
            selectCell  = laminaCells.findAt(leftBottomCorner)
            laminaPart.PartitionCellByPlanePointNormal(cells=
                selectCell,normal=laminaPart.datums[3].axis2,point=((-L/2)*0.99,(L/2)-((i+1)*d)-t,w))
    else:
        for i in range(N):
            laminaCells = laminaPart.cells
            selectCell  = laminaCells.findAt(rightTopCorner)
            laminaPart.PartitionCellByPlanePointNormal(cells=
                selectCell,normal=laminaPart.datums[3].axis2,point=((L/2),(-L/2)+((i+1)*d),w))
            laminaCells = laminaPart.cells
            selectCell  = laminaCells.findAt(rightTopCorner)
            laminaPart.PartitionCellByPlanePointNormal(cells=
                selectCell,normal=laminaPart.datums[3].axis2,point=((L/2),(-L/2)+((i+1)*d)+t,w))
#-----------------------------------------------------------------------------------------------
# Assigning Sections

for j in range(numLayers):
    laminaPart = laminaModel.parts['Ply-'+str(j+1)]
    if (Seq[j] > 0): 
        # a) Thin Strips
        for i in range(N+1):
           strip_point  = (((-L/2)+(i*d)+0.5*t),(-L/2)*0.99,w)
           laminaCells  = laminaPart.cells
           selectCell   = laminaCells.findAt(strip_point)
           strip_region = (selectCell,)
           laminaPart.SectionAssignment(region=strip_region, sectionName='MatrixSection')
        
        for i in range(N):
           strip_point  = ((-L/2)*0.99,((-L/2)+((i+1)*d)+(0.5*t)),w)
           laminaCells  = laminaPart.cells
           selectCell   = laminaCells.findAt(strip_point)
           strip_region = (selectCell,)
           laminaPart.SectionAssignment(region=strip_region, sectionName='MatrixSection')
        
        # b) Bulk Regions
        for i in range(N+1):
           bulk_point  = (((-L/2)+(i*d)+1.5*t),(-L/2)*0.99,w)
           laminaCells  = laminaPart.cells
           selectCell   = laminaCells.findAt(bulk_point)
           bulk_region = (selectCell,)
           laminaPart.SectionAssignment(region=bulk_region, sectionName='LaminaSection')
        
        for i in range(N+1):
           bulk_point  = ((-L/2)*0.99,((-L/2)+(i*d)+(1.5*t)),w)
           # bulk_point  = (L/2,((-L/2)+(i*d)+(1.5*t)),w)
           laminaCells  = laminaPart.cells
           selectCell   = laminaCells.findAt(bulk_point)
           bulk_region = (selectCell,)
           laminaPart.SectionAssignment(region=bulk_region, sectionName='LaminaSection')

    elif (Seq[j] < 0): 
        # a) Thin Strips
        for i in range(N+1):
           strip_point  = ( (L/2)*0.99,((-L/2)+(i*d)+0.5*t),w)
           laminaCells  = laminaPart.cells
           selectCell   = laminaCells.findAt(strip_point)
           strip_region = (selectCell,)
           laminaPart.SectionAssignment(region=strip_region, sectionName='MatrixSection')
        
        for i in range(N):
           strip_point  = ((-L/2)*0.99,((-L/2)+((i+1)*d)-(0.5*t)),w)
           laminaCells  = laminaPart.cells
           selectCell   = laminaCells.findAt(strip_point)
           strip_region = (selectCell,)
           laminaPart.SectionAssignment(region=strip_region, sectionName='MatrixSection')
        
        # b) Bulk Regions
        for i in range(N+1):
           bulk_point  = (((-L/2)+(i*d)+1.5*t),(-L/2)*0.99,w)
           laminaCells  = laminaPart.cells
           selectCell   = laminaCells.findAt(bulk_point)
           bulk_region = (selectCell,)
           laminaPart.SectionAssignment(region=bulk_region, sectionName='LaminaSection')
        
        for i in range(N+1):
           # bulk_point  = (-L/2,((-L/2)+(i*d)+(1.5*t)),w)
           bulk_point  = ((L/2)*0.99,((-L/2)+(i*d)+(1.5*t)),w)
           laminaCells  = laminaPart.cells
           selectCell   = laminaCells.findAt(bulk_point)
           bulk_region = (selectCell,)
           laminaPart.SectionAssignment(region=bulk_region, sectionName='LaminaSection')
    else:
        # a) Thin Strips
        for i in range(N):
            strip_point=((L/2),(-L/2)+((i+1)*d)+0.5*t,w)
            laminaCells  = laminaPart.cells
            selectCell   = laminaCells.findAt(strip_point)
            strip_region = (selectCell,)
            laminaPart.SectionAssignment(region=strip_region, sectionName='MatrixSection')
        for i in range(N+1):
            bulk_point  = ((L/2),((-L/2)+(i*d)+1.5*t),w)
            laminaCells  = laminaPart.cells
            selectCell   = laminaCells.findAt(bulk_point)
            bulk_region = (selectCell,)
            laminaPart.SectionAssignment(region=bulk_region, sectionName='LaminaSection')

#-----------------------------------------------------------------------------------------------
# Defining Interactions [Constraints]

# interface_pt1 = (0,0,w)
# MasterSurface1 = laminaModel.rootAssembly.instances['Ply-1-inst'].faces.findAt((interface_pt1,))
# laminaModel.rootAssembly.Surface(name='m_Surf-1', side1Faces=MasterSurface1)
# SlaveSurface1 = laminaModel.rootAssembly.instances['Ply-2-inst'].faces.findAt((interface_pt1,))
# laminaModel.rootAssembly.Surface(name='s_Surf-1', side1Faces=SlaveSurface1)
# laminaModel.Tie(adjust=ON, master= laminaModel.rootAssembly.surfaces['m_Surf-1'], name=
#     'Constraint-12', positionToleranceMethod=COMPUTED, slave=
#     laminaModel.rootAssembly.surfaces['s_Surf-1'], thickness=ON, tieRotations=ON)
# 
# interface_pt2 = (0,0,2*w)
# MasterSurface2 = laminaModel.rootAssembly.instances['Ply-2-inst'].faces.findAt((interface_pt2,))
# laminaModel.rootAssembly.Surface(name='m_Surf-2', side1Faces=MasterSurface2)
# SlaveSurface2 = laminaModel.rootAssembly.instances['Ply-3-inst'].faces.findAt((interface_pt2,))
# laminaModel.rootAssembly.Surface(name='s_Surf-2', side1Faces=SlaveSurface2)
# laminaModel.Tie(adjust=ON, master= laminaModel.rootAssembly.surfaces['m_Surf-2'], name=
#     'Constraint-23', positionToleranceMethod=COMPUTED, slave=
#     laminaModel.rootAssembly.surfaces['s_Surf-2'], thickness=ON, tieRotations=ON)
# 
# interface_pt3 = (0,0,3*w)
# MasterSurface3 = laminaModel.rootAssembly.instances['Ply-3-inst'].faces.findAt((interface_pt3,))
# laminaModel.rootAssembly.Surface(name='m_Surf-3', side1Faces=MasterSurface3)
# SlaveSurface3 = laminaModel.rootAssembly.instances['Ply-4-inst'].faces.findAt((interface_pt3,))
# laminaModel.rootAssembly.Surface(name='s_Surf-3', side1Faces=SlaveSurface3)
# laminaModel.Tie(adjust=ON, master= laminaModel.rootAssembly.surfaces['m_Surf-3'], name=
#     'Constraint-34', positionToleranceMethod=COMPUTED, slave=
#     laminaModel.rootAssembly.surfaces['s_Surf-3'], thickness=ON, tieRotations=ON)

#-----------------------------------------------------------------------------------------------
# Generating MESH

# Assign element type
# Linear element
elemType1 = mesh.ElemType(elemCode=C3D20R, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=C3D15, elemLibrary=STANDARD)
elemType3 = mesh.ElemType(elemCode=C3D10, elemLibrary=STANDARD)

for i in range(numLayers):
    p = laminaModel.parts['Ply-'+str(i+1)]
    c = p.cells
    f = p.faces
    mask_id = c.getMask() 
    pickedControl = c.getSequenceFromMask(mask=mask_id)
    p.Set(cells=pickedControl,name='Set-'+str(i+1))
    p.setMeshControls(elemShape=TET, regions=pickedControl, technique=FREE)
    p.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=meshSize)
    p.generateMesh()
    p.setElementType(regions=(pickedControl,) , elemTypes=(elemType1, elemType2, 
     elemType3))
    
#-----------------------------------------------------------------------------------------------
# Generate Cohesive Layers
# .
# .
# .

#-----------------------------------------------------------------------------------------------
# Create the Step
import step 

#Create a static general step
laminaModel.StaticStep(name='Load Step', previous='Initial',
   description='Apply boundary conditions in this step',
   nlgeom=ON)

#-----------------------------------------------------------------------------------------------
# Create the Field Output Request

laminaModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E'))

#-----------------------------------------------------------------------------------------------
# Create the History Output Request

# 'H-Output-Ply-i': History output of 'i'th' ply
# for j in range(numLayers):
#     laminaModel.HistoryOutputRequest(createStepName='Load Step', name=
#         'H-Output-Ply-'+str(j+1), rebar=EXCLUDE, region=
#         mdb.models['Lamina1'].rootAssembly.allInstances['Ply-'+str(j+1)+'-inst'].sets['Set-'+str(j+1)], 
#         sectionPoints=DEFAULT, variables=('S11', 'S22', 'S33', 'S12', 'S13', 
#             'S23', 'E11', 'E22', 'E33', 'E12', 'E13', 'E23'))

# laminaModel.HistoryOutputRequest(createStepName='Load Step', name=
#     'H-output-BD', rebar=EXCLUDE, region=laminaModel.rootAssembly.sets['BoundSet_CDGF'],
#     sectionPoints=DEFAULT, variables=('RF1','RF2'))

#-----------------------------------------------------------------------------------------------
# Boundary Conditions

BC_CDGF = laminaModel.DisplacementBC(amplitude=UNSET, createStepName=
    'Load Step', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=
    None, name='BC-CDGF', region=laminaModel.rootAssembly.sets['BoundSet_CDGF'], 
    u1=1.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

BC_B = laminaModel.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-B', 
    region=laminaModel.rootAssembly.sets['BoundSet_B'], u1=UNSET, u2=UNSET
    , u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

BC_A = laminaModel.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-A', 
    region=laminaModel.rootAssembly.sets['BoundSet_A'], u1=SET, u2=SET
    , u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

BC_BAHE = laminaModel.DisplacementBC(amplitude=UNSET, createStepName=
    'Load Step', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=
    None, name='BC-BAHE', region=laminaModel.rootAssembly.sets['BoundSet_BAHE'], 
    u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)

#-----------------------------------------------------------------------------------------------
# Creating Job
# NOTE: Update the subroutine path if applicable

# 'Job-with': with subroutine
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Lamina1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Job-with', nodalOutputPrecision=SINGLE, 
    numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
    ANALYSIS, userSubroutine=
    'C:\\ABAQUS+UMAT\\Week-31\\meka_orthoSCA3D.for'
    , waitHours=0, waitMinutes=0)

#-----------------------------------------------------------------------------------------------
# Saving the CAE file
mdb.saveAs(pathName='C:/ABAQUS+UMAT/Week-31/trial3_ver1.cae')





