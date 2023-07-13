from fenics import *
from classes.boundaries import *
from classes.createObjectSettings import *
from classes.problems import *
from classes.staticMethods import StaticMethod
from classes.thermalMaterial import ThermalParametersTransient, HeatConductivity
from classes.simulationParameters import SimulationParameters

#Testing steadystate solver
boxSpecimen = StaticMethod.createClassObjectFromJSON('settings/specimenSize.json', SpecimenSettings)
mesh = StaticMethod.createBoxMesh(boxSpecimen)
V = FunctionSpace(mesh, 'CG', 1)
temperature = StaticMethod.createClassObjectFromJSON('settings/temperatures.json', CreateObjectSettings)
top = Top(boxSpecimen)
bottom = Bottom()
bottomBoundaryTemperature = Constant(temperature.bottom)
topBoundaryTemperature = Constant(temperature.top)
temperatureField_n = interpolate(Constant(temperature.initial), V)
boundaryConditions = [DirichletBC(V, bottomBoundaryTemperature, bottom), 
                      DirichletBC(V, topBoundaryTemperature, top)]

#heatConductivity = HeatConductivity(temperatureDependency=True, temperatureDependentExpression="1+T", store=True)
heatConductivity = HeatConductivity(constant=100.0)
thermalParameters = ThermalParametersTransient(heatConductivity=heatConductivity)
simulationParameters = SimulationParameters()
transientHeatProblem = TransientHeat(V, initialTemperatureField=temperatureField_n, thermalParameters=thermalParameters, simulationParameters=simulationParameters)
transientHeatProblem.getSolution(boundaryConditions, nonZeroInitialGuess=True)
