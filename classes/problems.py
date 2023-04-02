from fenics import *
from abc import abstractstaticmethod, ABCMeta
from .outputFile import OutputFile
from .newFunction import NewFunction
from .thermalMaterial import ThermalParametersSteadyState, ThermalParametersTransient
from .simulationParameters import SimulationParameters

class IProblem(metaclass = ABCMeta):
    """Interface for thermal problems"""

    @abstractstaticmethod
    def defineVariationalProblem(self):
        """Interface function to define linear variational problem"""
        
    @abstractstaticmethod
    def getSolution(self):
        """Interface for solution to linear or non-linear problem"""


class SteadyStateHeat(IProblem):
    """Methods for solving steady state thermal problem"""
    def __init__(self, functionSpace, fileName="Steady-state", thermalParameters=None) -> None:
        self.functionSpace = functionSpace
        self.fileName = fileName
        if(thermalParameters==None):
            print("Assuming default thermal parameters")
            thermalParameters = ThermalParametersSteadyState()

        self.thermalParameters = thermalParameters
        if(thermalParameters.parameterSetType == "Linear-problem set"):
            self.solverType = "Linear"
        else:
            self.solverType = "Non-linear"

    def defineVariationalProblem(self, solution=None):
        v = TestFunction(self.functionSpace)
        du = TrialFunction(self.functionSpace)

        if(self.thermalParameters.heatConductivity.temperatureDependency == True):
            heatConductivity = self.thermalParameters.heatConductivityObject(solution)  #Returns a heat conductivity function
        else:
            heatConductivity = self.thermalParameters.heatConductivityObject
        
        if(self.thermalParameters.volumetricHeatSource.temperatureDependency == True):
            volumetricHeatSource = self.thermalParameters.volumetricHeatSourceObject(solution) #Returns a volumetric heat source function
        else:
            volumetricHeatSource = self.thermalParameters.volumetricHeatSourceObject

        if(self.solverType == 'Linear'):
            biLinearForm = inner(nabla_grad(v), heatConductivity*nabla_grad(du))*dx
            linearForm = volumetricHeatSource*v*dx
            
            return biLinearForm, linearForm
        
        elif(self.solverType == 'Non-linear'):
            if(solution == None):
                raise RuntimeError("Please provide the solution field variable as argument to this function")
            
            biLinearForm = inner(nabla_grad(v), heatConductivity*nabla_grad(solution))*dx
            linearForm = volumetricHeatSource*v*dx
            F = biLinearForm - linearForm
            J = derivative(F, solution, du)

            return F, J

    
    def solveSteadyStateProblem(self, boundaryConditions, solution, nonZeroInitialGuess):
        if (self.solverType == 'Linear'):
            print("Linear problem")
            biLinearForm, linearForm = self.defineVariationalProblem()
            solve(biLinearForm==linearForm, solution, boundaryConditions)
        
        elif self.solverType == 'Non-linear':
            print("Non-linear problem")
            F, J = self.defineVariationalProblem(solution)

            if(nonZeroInitialGuess == True):
                solution = self.solveWithNonZeroInitialization(solution)

            problem = NonlinearVariationalProblem(F, solution, boundaryConditions, J)
            solver = NonlinearVariationalSolver(problem)
            solver.solve()


    def getSolution(self, boundaryConditions, nonZeroInitialGuess=False):
        """Combined function to solve linear and non-linear FEM problem"""

        solution = NewFunction(self.functionSpace, name="Temperature", store=True)
        listOfFields = {"Solution":solution}   #Default list of fields
        listOfFields = self.setupListOfStoredFields(listOfFields)

        file = OutputFile(self.fileName, listOfFields)

        if(self.solverType == 'Linear' or self.solverType == 'Non-linear'):
            self.solveSteadyStateProblem(boundaryConditions, solution, nonZeroInitialGuess)
        
        else:
            raise Exception("solverType can be Linear or Non-linear")

        self.assignListOfStoredFields(listOfFields)
        file.exportFile() #Exporting fields to output file



    def assignListOfStoredFields(self, listOfFields):

        if(self.thermalParameters.heatConductivity.store == True or self.thermalParameters.volumetricHeatSource.store == True):
            if(self.thermalParameters.heatConductivity.temperatureDependentTable != None or 
               self.thermalParameters.volumetricHeatSource.temperatureDependentTable != None):
                bufferSolution = NewFunction(self.functionSpace)
                bufferSolution.assign(listOfFields["Solution"])
                
        if(self.thermalParameters.heatConductivity.store == True):
            if(self.thermalParameters.heatConductivity.temperatureDependentTable != None):
                heatConductivity = self.thermalParameters.heatConductivityObject(bufferSolution)
                listOfFields["Heat Conductivity"].assign(heatConductivity)
            elif(self.thermalParameters.heatConductivity.temperatureDependentExpression != None):
                heatConductivityVector = self.thermalParameters.heatConductivityObject(listOfFields["Solution"].vector()[:])
                listOfFields["Heat Conductivity"].vector()[:] = heatConductivityVector
            else:
                heatConductivity = interpolate(self.thermalParameters.heatConductivityObject, self.functionSpace)
                listOfFields["Heat Conductivity"].assign(heatConductivity)

        if(self.thermalParameters.volumetricHeatSource.store == True):
            if(self.thermalParameters.volumetricHeatSource.temperatureDependentTable != None):
                volumetricHeatSourceFieldBuffer = self.thermalParameters.volumetricHeatSourceObject(bufferSolution)
                listOfFields["Volumetric Heat Source"].assign(volumetricHeatSourceFieldBuffer)
            elif(self.thermalParameters.volumetricHeatSource.temperatureDependentExpression != None):
                volumetricHeatSourceVector = self.thermalParameters.volumetricHeatSourceObject(listOfFields["Solution"].vector()[:])
                listOfFields["Volumetric Heat Source"].vector()[:] = volumetricHeatSourceVector
            else:
                volumetricHeatSourceField = interpolate(self.thermalParameters.volumetricHeatSourceObject, self.functionSpace)
                listOfFields["Volumetric Heat Source"].assign(volumetricHeatSourceField)


    def setupListOfStoredFields(self, listOfFields):

        if self.thermalParameters.heatConductivity.store == True:
            heatConductivityField = NewFunction(self.functionSpace, name="Heat Conductivity")
            heatConductivityField.store = True
            listOfFields["Heat Conductivity"] = heatConductivityField

        if self.thermalParameters.volumetricHeatSource.store == True:
            volumetricHeatSourceField = NewFunction(self.functionSpace, name="Volumetric Heat Source")
            volumetricHeatSourceField.store = True
            listOfFields["Volumetric Heat Source"] = volumetricHeatSourceField

        return listOfFields

    def solveWithNonZeroInitialization(self, solution):
        firstGuess = 1.0     #Non zero value 
        solution.vector()[:] = firstGuess
        parameters["krylov_solver"]["nonzero_initial_guess"] = True
        return solution
    

class TransientHeat(SteadyStateHeat):
    def __init__(self, functionSpace, initialTemperatureField, fileName="transient", thermalParameters=None, simulationParameters=None) -> None:
        if thermalParameters == None:
            thermalParameters = ThermalParametersTransient()
            
        super().__init__(functionSpace, fileName, thermalParameters)

        if simulationParameters is None:
            simulationParameters = SimulationParameters()

        self.simulationParameters = simulationParameters
        self.initialTemperatureField = initialTemperatureField


    def defineVariationalProblem(self, T_previous, solution=None):
        v = TestFunction(self.functionSpace)
        du = TrialFunction(self.functionSpace)
        theta = self.simulationParameters.trapezoidalParameter
        dt = self.simulationParameters.dt
        dynamicFieldStorage = False

        if(self.thermalParameters.heatConductivity.temperatureDependency == True):
            kappa = self.thermalParameters.heatConductivityObject(solution)  #Returns a heat conductivity function
            dynamicFieldStorage = True
        else:
            kappa = self.thermalParameters.heatConductivityObject
        
        if(self.thermalParameters.volumetricHeatSource.temperatureDependency == True):
            f = self.thermalParameters.volumetricHeatSourceObject(solution) #Returns a volumetric heat source function
            dynamicFieldStorage = True
        else:
            f = self.thermalParameters.volumetricHeatSourceObject

        if(self.thermalParameters.density.temperatureDependency == True):
            rho = self.thermalParameters.densityObject(solution) #Returns a density function
            dynamicFieldStorage = True
        else:
            rho = self.thermalParameters.densityObject

        if(self.thermalParameters.heatCapacity.temperatureDependency == True):
            c = self.thermalParameters.heatCapacityObject(solution) #Returns a density function
            dynamicFieldStorage = True
        else:
            c = self.thermalParameters.heatCapacityObject

        if(self.solverType == 'Linear'):
            biLinearForm = theta*dt*inner(nabla_grad(v), kappa*nabla_grad(du))*dx + rho*c*du*v*dx
            linearForm = (rho*c*T_previous*v + dt*f*v - (1-theta)*dt*inner(nabla_grad(v), kappa*nabla_grad(du)))*dx

            return biLinearForm, linearForm

        elif(self.solverType == 'Non-linear'):
            if(solution == None):
                raise RuntimeError("Pass the solution field as argument in this function")

            biLinearForm = (theta*dt*inner(nabla_grad(v), kappa*nabla_grad(solution)) + rho*c*solution*v)*dx
            linearForm = (rho*c*T_previous*v + dt*f*v - (1-theta)*dt*inner(nabla_grad(v), kappa*nabla_grad(solution)))*dx

            F = biLinearForm - linearForm
            J = derivative(F, solution, du)

            return F, J, dynamicFieldStorage
            

    def getSolution(self, boundaryConditions, nonZeroInitialGuess=False):
        """Combined function to solve linear and non-linear FEM problem"""

        T_previous = self.initialTemperatureField
        solution = NewFunction(self.functionSpace, name="Temperature", store=True)
        listOfFields = {"Solution": solution}   #Default list of fields

        listOfFields = self.setupListOfStoredFields(listOfFields)   #Add all the fields to be stored to the listOfFields
        file = OutputFile(self.fileName, listOfFields)

        if(self.solverType == 'Linear' or self.solverType =='Non-linear'):
            if(self.simulationParameters.trapezoidalParameter==1):
                print("Implicit time stepping ---- for transient thermal problem")
            elif(self.simulationParameters.trapezoidalParameter==0):
                print("Explicit time stepping ---- for transient thermal problem")
            else:
                print("Crank-Nicolson time stepping (alpha =", str(self.simulationParameters.trapezoidalParameter)+")",
                        "---- for transient thermal problem")

            self.solveTransientProblem(T_previous, solution, listOfFields, file, boundaryConditions, nonZeroInitialGuess)
        
        else:
            raise Exception("solverType can be Linear or Non-linear")


    def setupListOfStoredFields(self, listOfFields):
        listOfFields = super().setupListOfStoredFields(listOfFields)

        if self.thermalParameters.density.store == True:
            densityField = NewFunction(self.functionSpace, name="Density")
            densityField.store = True
            listOfFields["Density"] = densityField

        if self.thermalParameters.heatCapacity.store == True:
            heatCapacityField = NewFunction(self.functionSpace, name="Heat Capacity")
            heatCapacityField.store = True
            listOfFields["Heat Capacity"] = heatCapacityField

        return listOfFields
    
    def assignListOfStoredFields(self, listOfFields):
        super().assignListOfStoredFields(listOfFields)

        if(self.thermalParameters.density.store is True or self.thermalParameters.heatCapacity.store is True):
            if(self.thermalParameters.density.temperatureDependentTable != None or 
               self.thermalParameters.heatCapacity.temperatureDependentTable != None):
                bufferSolution = NewFunction(self.functionSpace)
                bufferSolution.assign(listOfFields["Solution"])
                
        if(self.thermalParameters.density.store is True):
            if(self.thermalParameters.density.temperatureDependentTable != None):
                density = self.thermalParameters.densityObject(bufferSolution)
                listOfFields["Density"].assign(density)
            elif(self.thermalParameters.density.temperatureDependentExpression != None):
                densityVector = self.thermalParameters.densityObject(listOfFields["Solution"].vector()[:])
                listOfFields["Density"].vector()[:] = densityVector
            else:
                density = interpolate(self.thermalParameters.densityObject, self.functionSpace)
                listOfFields["Density"].assign(density)

        if(self.thermalParameters.heatCapacity.store is True):
            if(self.thermalParameters.heatCapacity.temperatureDependentTable != None):
                heatCapacity = self.thermalParameters.heatCapacityObject(bufferSolution)
                listOfFields["Heat Capacity"].assign(heatCapacity)
            elif(self.thermalParameters.heatCapacity.temperatureDependentExpression != None):
                heatCapacityVector = self.thermalParameters.heatCapacityObject(listOfFields["Solution"].vector()[:])
                listOfFields["Heat Capacity"].vector()[:] = heatCapacityVector
            else:
                heatCapacity = interpolate(self.thermalParameters.heatCapacityObject, self.functionSpace)
                listOfFields["Heat Capacity"].assign(heatCapacity)


    def solveTransientProblem(self, T_previous, solution, listOfFields, file, boundaryConditions, nonZeroInitialGuess):
        totalTime =  self.simulationParameters.totalTime
        dt =  self.simulationParameters.dt
        timeTolerance = self.simulationParameters.timeTolerance
        storeFrequency = self.simulationParameters.storeFrequency
        timeStepCount = 0
        dynamicFieldStorage = False
        
        if self.solverType == 'Linear':
            biLinearForm, linearForm = self.defineVariationalProblem(T_previous)  
            
        elif self.solverType == 'Non-linear':
            F, J, dynamicFieldStorage = self.defineVariationalProblem(T_previous, solution)

            if(dynamicFieldStorage == False):
                self.assignListOfStoredFields(listOfFields)


        currentTime = 0
        storedFrameCount = 0

        while currentTime < totalTime:
            timeStepCount += 1
            currentTime += dt
            if((abs(currentTime-totalTime)) < timeTolerance):
                currentTime = totalTime

            storeFrame = (timeStepCount%storeFrequency == 0) or (currentTime == totalTime)
            if storeFrame:
                storedFrameCount += 1
                if(storedFrameCount == 1):
                    print("Storing frames at:")
                set_log_active(True)
                print("Time :", str(round(currentTime, 6)), "  Time-Step :", str(timeStepCount))
            else:
                set_log_active(False)

            if self.solverType == 'Linear':
                solve(biLinearForm==linearForm, solution, boundaryConditions)    
            
            elif self.solverType == 'Non-linear':
                if(nonZeroInitialGuess == True):
                    solution = self.solveWithNonZeroInitialization(solution)

                problem = NonlinearVariationalProblem(F, solution, boundaryConditions, J)
                solver = NonlinearVariationalSolver(problem)
                solver.solve()
                if(dynamicFieldStorage == True):
                    self.assignListOfStoredFields(listOfFields)    

            else: 
                raise Exception("solverType can be Linear or Non-linear")
            
            if storeFrame:
                file.exportFile(currentTime)       #Exporting fields to output file

            T_previous.assign(solution)


class TransientMovingHeatSource(IProblem):
    """To be implemented in a clean way by inheriting form steadystate class"""
    
    def __init__(self, thermalMaterialParameters, printSettings, fileName) -> None:
        self.density = thermalMaterialParameters.density
        self.heatCapacity = thermalMaterialParameters.heatCapacity
        self.heatConductivity = thermalMaterialParameters.heatConductivity
        self.fileName = fileName
        self.specimenLength = printSettings.length
        self.specimenWidth = printSettings.width
        self.specimenThickness = printSettings.thickness
        self.hsX = printSettings.heatSourceSizeX / 1000   #in mmm
        self.hsY = printSettings.heatSourceSizeY / 1000   #in mmm
        self.hsZ = printSettings.heatSourceSizeZ / 1000   #in mmm
        self.scanSpeed = printSettings.scanSpeed
        self.totalTime = printSettings.length / printSettings.scanSpeed
        self.heatInput = (printSettings.laserEfficiency*printSettings.laserPower)/(self.hsX*self.hsY*self.hsZ) 
        self.dt = 1e-5 
        self.storeFrequency = 10

    def getCubicMovingHeatSource(self, time):
        heatSource = Expression("(x[0] > v*t & x[0] < hsX + v*t) & (x[1] > width/2 - hsZ/2 & x[1] < width/2 + hsZ/2) & (x[2] > thickness - hsY) ?  q : 0", 
                                v = self.scanSpeed, t = time, hsX = self.hsX, thickness = self.specimenThickness, hsY = self.hsY, 
                                width = self.specimenWidth, hsZ = self.hsZ, q = self.heatInput, degree = 0)
        return heatSource
        
    def defineVariationalProblem(self, functionSpace, previousTemperature, heatSource):
        self.v = TestFunction(functionSpace)
        self.u = TrialFunction(functionSpace)
        self.F = self.density*self.heatCapacity*self.u*self.v*dx + self.dt*self.heatConductivity*dot(grad(self.u), grad( self.v))*dx - (self.density*self.heatCapacity*previousTemperature + self.dt*heatSource)* self.v*dx  
        self.a, self.L =  lhs(self.F), rhs(self.F)
        return self.a, self.L

    def linearSolve(self, functionSpace, boundaryConditions, previousTemperature):
        fileObject = File(self.fileName)
        solution = Function(functionSpace, name = "Temperature")
        solution.assign(previousTemperature)
        currentTime = 0.0
        fileObject<<(solution, currentTime)
        heatSource = self.getCubicMovingHeatSource(currentTime)
        timeStep = 0        
        biLinearForm, linearForm = self.defineVariationalProblem(functionSpace, previousTemperature, heatSource)

        while currentTime < self.totalTime:
            timeStep += 1
            solve(biLinearForm==linearForm, solution, boundaryConditions)
            if(timeStep%self.storeFrequency == 0):
                print("time (steps: " + str(timeStep) + " time:" + str(currentTime))
                fileObject<< (solution, currentTime)
            previousTemperature.assign(solution)
            currentTime += self.dt
            heatSource.t = currentTime


    def defineNonLinearVariationalProblem(self):
        raise RuntimeError("Not implemented")
    
    def nonLinearSolve(self):
        raise RuntimeError("Not implemented")


    def checkEnergyBalance(self):
        raise RuntimeError("Not implemented")

        
    def getSolution(self):
        raise RuntimeError("Not implemented")
 