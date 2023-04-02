import json
from fenics import *


class StaticMethod:
    @staticmethod
    def createClassObjectFromJSON(filename, className):
        with open(filename) as var:
            var2 = json.load(var)

        classInstance = className(**var2)
        return classInstance

    @staticmethod
    def createBoxMesh(specimenSetting):
        meshObject = BoxMesh(Point(0.0, 0.0, 0.0), Point(specimenSetting.length, specimenSetting.thickness, specimenSetting.width), 
                             specimenSetting.cellCountX, specimenSetting.cellCountY, specimenSetting.cellCountZ)
        return meshObject
    
    @staticmethod
    def createBoxMeshFromPrintSettings(printSetting):
        meshObject = BoxMesh(Point(0.0, 0.0, 0.0), Point(printSetting.length, printSetting.width, printSetting.thickness), printSetting.cellCountX, printSetting.cellCountY, printSetting.cellCountZ)
        return meshObject
    
    