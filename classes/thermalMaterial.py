from fenics import *
from array import *
import numpy as np
from .volumetricHeatSource import VolumetricHeatSource, MovingVolumetricHeatSource
from .heatConductivity import HeatConductivity
from .density import Density
from .heatCapacity import HeatCapacity


class ThermalParametersSteadyState:

    """Class containing information and methods for thermal parameters"""
 
    def __init__(self, heatConductivity=None, volumetricHeatSource=None) -> None:

        if(heatConductivity == None):
            print("Heat Conductivity not provided.")
            print("Assuming constant Heat conductivity = 1.0.")
            heatConductivity = HeatConductivity()
              
        if(volumetricHeatSource == None):
            print("Heat Source not provided.")
            print("Assuming constant Heat Source = 0.0.")
            volumetricHeatSource = VolumetricHeatSource()

        self.heatConductivity = heatConductivity
        self.volumetricHeatSource = volumetricHeatSource
        self.parameterSetType = "Linear-problem set"
        
        self.heatConductivityObject = Constant(self.heatConductivity.constant)
        self.volumetricHeatSourceObject = Constant(self.volumetricHeatSource.constant)

        if(heatConductivity.spaceDependency == True):
            self.heatConductivityObject = Expression(heatConductivity.spatialExpression, degree=0)
            self.parameterSetType = "Non-linear-problem set"

        if(volumetricHeatSource.spaceDependency == True):
            self.volumetricHeatSourceObject = Expression(volumetricHeatSource.spatialExpression, degree=0)
            if(heatConductivity.spaceDependency == True or heatConductivity.temperatureDependency == True):
                self.parameterSetType = "Non-linear-problem set"
            else:
                self.parameterSetType = "Linear-problem set"
        
        if(heatConductivity.temperatureDependency == True or volumetricHeatSource.temperatureDependency == True):
            self.parameterSetType = "Non-linear-problem set"

            if(heatConductivity.temperatureDependency == True):
                if(heatConductivity.temperatureDependentExpression != None):
                    self.heatConductivityObject = self.heatConductivityFunction
                elif(heatConductivity.temperatureDependentTable != None):
                    self.heatConductivityObject = self.heatConductivityFromTable
                
            if(volumetricHeatSource.temperatureDependency == True):
                if(volumetricHeatSource.temperatureDependentExpression != None):
                    self.volumetricHeatSourceObject = self.volumetricHeatSourceFunction
                elif(volumetricHeatSource.temperatureDependentTable != None):
                    self.volumetricHeatSourceObject = self.volumetricHeatSourceFromTable

                
    def heatConductivityFunction(self, variable):
        """Takes string of mathematical expression in terms of temperature T, returns the mathematical expression"""
        return eval(self.heatConductivity.temperatureDependentExpression, {"T": variable})
    
    def volumetricHeatSourceFunction(self, variable):
        """Takes string of mathematical expression in terms of temperature T, returns the mathematical expression"""
        return eval(self.volumetricHeatSource.temperatureDependentExpression, {"T": variable})
    
    def heatConductivityFromTable(self, variable):
        """Input table as a 2D array from CSV file, output: linearly interpolated value"""
    
        CSVData = open(self.heatConductivity.temperatureDependentTable)  #this is a csv file
        variable = self.getInterpolatedVector(variable, CSVData)

        return variable

    
    def volumetricHeatSourceFromTable(self, variable):
        """Input table as a 2D array from CSV file, output: linearly interpolated value"""
        CSVData = open(self.volumetricHeatSource.temperatureDependentTable)  #this is a csv file
        variable = self.getInterpolatedVector(variable, CSVData)

        return variable
    
    def getInterpolatedVector(self, variable, CSVData):
        temperatureParameterArray = np.genfromtxt(CSVData, delimiter=",")
        variableVector = variable.vector()[:]
        parameterVector = len(variableVector)*[0]
        i = 0
        for eachVectorEntry in variableVector:
            parameterVector[i] = self.getInterpolatedValue(eachVectorEntry, temperatureParameterArray)
            i += 1

        variable.vector()[:] = parameterVector[:]

        return variable


        
    def getInterpolatedValue(self, temperature, temperatureParameterArray):
        """To be used in the function heatConductivityFromTable, not in use currently. """

        if(len(temperatureParameterArray) < 2):
            raise RuntimeError("Incorrect data table for material parameter.")

        for i in range(len(temperatureParameterArray)):
            if(i == 0):
                if(temperature < temperatureParameterArray[i][0]):
                    conductivityValue = temperatureParameterArray[i][1]
                else:
                    conductivityValue = self.interpolate(temperature, temperatureParameterArray[i][0], temperatureParameterArray[i+1][0], 
                                                    temperatureParameterArray[i][1], temperatureParameterArray[i+1][1])

            elif(i == len(temperatureParameterArray) - 1):
                if(temperature > temperatureParameterArray[i][0]):
                    conductivityValue = temperatureParameterArray[i][1]
                else:
                    conductivityValue = self.interpolate(temperature, temperatureParameterArray[i-1][0], temperatureParameterArray[i][0], 
                                                    temperatureParameterArray[i-1][1], temperatureParameterArray[i][1])

            elif(temperature >= temperatureParameterArray[i][0] and temperature <= temperatureParameterArray[i+1][0]):
                conductivityValue = self.interpolate(temperature, temperatureParameterArray[i][0], temperatureParameterArray[i+1][0], 
                                                    temperatureParameterArray[i][1], temperatureParameterArray[i+1][1])

            else:
                print("Problem with index", str(i))
                raise RuntimeError("There is something wrong in this function")

        return conductivityValue
    
    
    def interpolate(self, temperature, temperature1, temperature2, parameter1, parameter2):
        slope = (parameter2 - parameter1)/(temperature2 - temperature1)
        interpolatedValue = parameter1 + slope*(temperature - temperature1)
        return interpolatedValue



class ThermalParametersTransient(ThermalParametersSteadyState):
    def __init__(self, density=None, heatCapacity=None, heatConductivity=None, volumetricHeatSource=None) -> None:
        super().__init__(heatConductivity, volumetricHeatSource)

        if(density is None):
            print("Density not provided.")
            print("Assuming constant density = 1.0.")
            density = Density()

        if(heatCapacity is None):
            print("Heat capacity not provided.")
            print("Assuming constant Heat capacity = 1.0.")
            heatCapacity = HeatCapacity()

        self.density = density
        self.heatCapacity = heatCapacity

        self.densityObject = Constant(self.density.constant)
        self.heatCapacityObject = Constant(self.heatCapacity.constant)

        if(density.spaceDependency == True or heatCapacity.spaceDependency == True):
            self.parameterSetType = "Non-linear-problem set"

            if(density.spaceDependency == True):
                self.densityObject = Expression(density.spatialExpression, degree=0)

            if(heatCapacity.spaceDependency == True):
                self.heatCapacityObject = Expression(heatCapacity.spatialExpression, degree=0)


        if(density.temperatureDependency == True or heatCapacity.temperatureDependency == True):
            self.parameterSetType = "Non-linear-problem set"

            if (density.temperatureDependency == True):
                if(density.temperatureDependentExpression != None):
                    self.densityObject = self.densityFunction
                elif(density.temperatureDependentTable != None):
                    self.densityObject = self.densityFromTable

            if (heatCapacity.temperatureDependency == True):
                if(heatCapacity.temperatureDependentExpression != None):
                    self.heatCapacityObject = self.heatCapacityFunction
                elif(heatCapacity.temperatureDependentTable != None):
                    self.heatCapacityObject = self.heatCapacityFromTable


    def densityFunction(self, variable):
        """Takes string of mathematical expression in terms of temperature T, returns the mathematical expression"""
        return eval(self.density.temperatureDependentExpression, {"T": variable})
    

    def densityFromTable(self, variable):
        """Input table as a 2D array from CSV file, output: linearly interpolated value"""
        CSVData = open(self.density.temperatureDependentTable)  #this is a csv file
        variable = self.getInterpolatedVector(variable, CSVData)

        return variable
    

    def heatCapacityFunction(self, variable):
        """Takes string of mathematical expression in terms of temperature T, returns the mathematical expression"""
        return eval(self.heatCapacity.temperatureDependentExpression, {"T": variable})
    

    def heatCapacityFromTable(self, variable):
        """Input table as a 2D array from CSV file, output: linearly interpolated value"""
        CSVData = open(self.heatCapacity.temperatureDependentTable)  #this is a csv file
        variable = self.getInterpolatedVector(variable, CSVData)

        return variable


class ThermalParametersTransientLaserBeam(ThermalParametersTransient):
    def __init__(self, printSettings, density=None, heatCapacity=None, heatConductivity=None, volumetricHeatSource=None) -> None:
        if volumetricHeatSource != None:
            print("Not using the input volumetric heat source")
            print("Computing moving volumetric heat source from print settings")

        volumetricHeatSource = MovingVolumetricHeatSource(printSettings)  #Reassign volumetric heat source
        super().__init__(density, heatCapacity, heatConductivity, volumetricHeatSource)

        self.volumetricHeatSourceObject = self.getCubicMovingHeatSource(volumetricHeatSource.printSettings)
         

    def getCubicMovingHeatSource(self, printSettings):
        heatSource = Expression("(x[0] > v*t & x[0] < hsX + v*t) & (x[1] > width/2 - hsZ/2 & x[1] < width/2 + hsZ/2) & (x[2] > thickness - hsY) ?  q : 0", 
                                v = printSettings.scanSpeed, t = 0.0, hsX = printSettings.cubeEdgeLength, thickness = printSettings.thickness, hsY = printSettings.cubeEdgeLength, 
                                width = printSettings.width, hsZ = printSettings.cubeEdgeLength, q = printSettings.laserHeatInput, degree = 0)
        return heatSource

