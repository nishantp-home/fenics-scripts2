from .printerSpecifications import PrinterSpecifications
from .specimenSize import BoxSpecimen
from fenics import BoxMesh, Point

class PrintSettings(BoxSpecimen, PrinterSpecifications):
    """Class containing data for print simulation"""
    
    cubeEdgeLength = 0.05               #Incorporate heat source size in the class volumetricHeatSource
    __cellsPerHeatSourceSide = 2
    __cellsPerUnitHeatSourceSide = __cellsPerHeatSourceSide/cubeEdgeLength

    
    def __init__(self, length=0.2, width=0.1, thickness=0.5, cellsPerUnitDimensions=10, 
                 scanSpeed=1000, laserPower=300, laserEfficiency=0.7, hatchDistance=0.05) -> None:
        BoxSpecimen.__init__(self, length, width, thickness, cellsPerUnitDimensions)
        PrinterSpecifications.__init__(self, scanSpeed, laserPower, laserEfficiency, hatchDistance)
        self.laserHeatInput = self.__getLaserHeatInput()
        self.__cellCountX = int(self.length*self.__cellsPerUnitHeatSourceSide)
        self.__cellCountY = int(self.thickness*self.__cellsPerUnitHeatSourceSide)
        self.__cellCountZ= int(self.width*self.__cellsPerUnitHeatSourceSide)


    def __getLaserHeatInput(self):
        heatInput = (self.laserPower*self.laserEfficiency)/(self.cubeEdgeLength**3)   #get HeatSource Size here
        return heatInput
    
    def createBoxMesh(self):
         meshObject = BoxMesh(Point(0.0, 0.0, 0.0), Point(self.length, self.thickness, self.width), 
                             self.__cellCountX, self.__cellCountY, self.__cellCountZ)
         
         return meshObject






        