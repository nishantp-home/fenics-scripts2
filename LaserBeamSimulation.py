from fenics import *
from classes.boundaries import *
from classes.createObjectSettings import *
from classes.problems import *
from classes.staticMethods import StaticMethod
from classes.thermalMaterial import ThermalParametersTransientLaserBeam, HeatCapacity, HeatConductivity
from classes.printSettings import PrintSettings


# Print settings to setup 3D print Geometry and printer specifications
printSettings = PrintSettings(length=1.0, width=0.3, thickness=0.3)     #Required

# Generate FEM mesh from print settings
# Generating a mesh out of a 3D box of a given dimension (i.e. length, width, thickness) using the createBoxMesh method
# Alternatively, a mesh could also be imported here
mesh = printSettings.createBoxMesh()    #Required

# Function space from mesh
V = FunctionSpace(mesh, 'CG', 1)    #Required


# Initial and boundary conditions for temperature 
# Currently input from JSON files
# This need to be modularized
temperature = StaticMethod.createClassObjectFromJSON('settings/temperatures.json', CreateObjectSettings)   #Required
bottom = Bottom()  #Required if bottom boundary is Dirichlet
 
bottomBoundaryTemperature = Constant(temperature.bottom)

# Initial temperature field over the mesh
temperatureField_n = interpolate(Constant(temperature.initial), V)    #Required for transient thermal problems

boundaryConditions = [DirichletBC(V, bottomBoundaryTemperature, bottom)]  #Required if there are Dirichlet boundaries



# Optional thermal parameters setup
# If not provided, they take default constant values provided in their respective class constructors
heatCapacity = HeatCapacity(constant=0.0033)
heatConductivity = HeatConductivity(constant=0.05)
thermalParameters = ThermalParametersTransientLaserBeam(printSettings=printSettings, heatCapacity=heatCapacity, heatConductivity=heatConductivity)


# Main problem setup and solution ###############################################
printProblem = TransientHeatLaserBeam(V, initialTemperatureField=temperatureField_n, printSettings=printSettings, thermalParameters=thermalParameters)

printProblem.getSolution(boundaryConditions, nonZeroInitialGuess=True)
#################################################################################