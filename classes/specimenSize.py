from fenics import *

class BoxSpecimen:

    """"Dimensions of box specimen"""
    def __init__(self, length=0.2, width=0.1, thickness=0.5, cellsPerUnitDimensions=10) -> None:
        self.length = length
        self.width = width
        self.thickness = thickness

        self.__cellsPerUnitDimensions = cellsPerUnitDimensions
        self.__cellCountX = int(self.length*self.__cellsPerUnitDimensions)
        self.__cellCountY = int(self.thickness*self.__cellsPerUnitDimensions)
        self.__cellCountZ = int(self.width*self.__cellsPerUnitDimensions)

    def createBoxMesh(self):
         meshObject = BoxMesh(Point(0.0, 0.0, 0.0), Point(self.length, self.thickness, self.width), 
                             self.__cellCountX, self.__cellCountY, self.__cellCountZ)
         
         return meshObject
