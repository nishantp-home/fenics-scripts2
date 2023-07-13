from fenics import *
from classes.boundaries import *
from classes.createObjectSettings import *
from classes.problems import *
from classes.staticMethods import StaticMethod
from classes.thermalMaterial import ThermalParametersSteadyState, VolumetricHeatSource, HeatConductivity
import numpy as np


#Testing steadystate solver
boxSpecimen = StaticMethod.createClassObjectFromJSON('settings/specimenSize.json', SpecimenSettings)
mesh = StaticMethod.createBoxMesh(boxSpecimen)
V = FunctionSpace(mesh, 'CG', 1)


bc = MYBC(boxSpecimen, mesh, bottom="D", top="D")








temperature = StaticMethod.createClassObjectFromJSON('settings/temperatures.json', CreateObjectSettings)
top = Top(boxSpecimen)
bottom = Bottom()
bottomBoundaryTemperature = Constant(temperature.bottom)
topBoundaryTemperature = Constant(temperature.top)

#################################################
boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
boundary_markers.set_all(7)
bottom.mark(boundary_markers, 0)
top.mark(boundary_markers, 1)
 


class BoundaryConditions:
    def __init__(self, bcType, marker, valueObject) -> None:
        self.bcType = bcType
        self.marker = marker
        self.valueObject = valueObject

class BcList:
    def __init__(self, mesh, BCObjectList) -> None:
        self.boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1).set_all(7)
        self.BCObjectList = BCObjectList

    def getDirichletBCs(self):
        boundaryConditions = []
        for currentBC in self.BCObjectList:
            if currentBC.bcType == "Dirichlet":
                bc = DirichletBC(V, currentBC.valueObject, self.boundary_markers, currentBC.marker)
                boundaryConditions.append(bc)

        return boundaryConditions



boundary_conditions = [BoundaryConditions("Dirichlet", 0, bottomBoundaryTemperature),
                       BoundaryConditions("Dirichlet", 1, topBoundaryTemperature)]


bcs = BcList(mesh, boundary_conditions)
boundaryConditions = bcs.getDirichletBCs()



# boundaryConditions = []
# for currentBC in boundary_conditions:
#     if currentBC.bcType == "Dirichlet":
#         bc = DirichletBC(V, currentBC.valueObject,
#                          boundary_markers, currentBC.marker)
#         boundaryConditions.append(bc)


# boundary_conditions = {0: {'Dirichlet': bottomBoundaryTemperature},
                    #    1: {'Dirichlet': topBoundaryTemperature}}

# boundaryConditions = []
# for i in boundary_conditions:
#     if 'Dirichlet' in boundary_conditions[i]:
#         bc = DirichletBC(V, boundary_conditions[i]['Dirichlet'],
#                          boundary_markers, i)
#         boundaryConditions.append(bc)





###################################################

heatConductivity = HeatConductivity()
volumetricHeatSource = VolumetricHeatSource()
thermalParameters = ThermalParametersSteadyState(heatConductivity=heatConductivity, volumetricHeatSource=volumetricHeatSource)

# steadyStateHeatProblem = SteadyStateHeat(V, thermalParameters=thermalParameters)
# steadyStateHeatProblem.getSolution(boundaryConditions)


ssproblem = SteadyStateNew(V, boundaryConditions=boundaryConditions, thermalParameters=thermalParameters)
ssproblem.getSolution()

