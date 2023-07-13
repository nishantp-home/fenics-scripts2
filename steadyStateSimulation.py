from fenics import *
from classes.boundaries import *
from classes.createObjectSettings import *
from classes.problems import *
from classes.staticMethods import StaticMethod
from classes.thermalMaterial import ThermalParametersSteadyState, VolumetricHeatSource, HeatConductivity

#Testing steadystate solver
boxSpecimen = StaticMethod.createClassObjectFromJSON('settings/specimenSize.json', SpecimenSettings)
mesh = StaticMethod.createBoxMesh(boxSpecimen)
V = FunctionSpace(mesh, 'CG', 1)
temperature = StaticMethod.createClassObjectFromJSON('settings/temperatures.json', CreateObjectSettings)
top = Top(boxSpecimen)
bottom = Bottom()
bottomBoundaryTemperature = Constant(temperature.bottom)
topBoundaryTemperature = Constant(temperature.top)
boundaryConditions = [DirichletBC(V, bottomBoundaryTemperature, bottom), 
                      DirichletBC(V, topBoundaryTemperature, top)]

heatConductivity = HeatConductivity(spaceDependency=True, spatialExpression="1+x[1]", store=True)
volumetricHeatSource = VolumetricHeatSource(1.0, store=True)
thermalParameters = ThermalParametersSteadyState(heatConductivity=heatConductivity, volumetricHeatSource=volumetricHeatSource)

steadyStateHeatProblem = SteadyStateHeat(V, thermalParameters=thermalParameters)
steadyStateHeatProblem.getSolution(boundaryConditions)
