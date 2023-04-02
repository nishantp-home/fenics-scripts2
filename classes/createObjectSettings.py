class CreateObjectSettings:
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            if type(v) is dict:
                setattr(self, k, CreateObjectSettings(**v))
            else:
                setattr(self, k, v)

class SimulationSettings(CreateObjectSettings):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.dt = self.totalTime/self.numTimeSteps

class SpecimenSettings(CreateObjectSettings):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.cellCountX = int(self.length*self.cellsPerUnitDimension)
        self.cellCountY = int(self.thickness*self.cellsPerUnitDimension)
        self.cellCountZ = int(self.width*self.cellsPerUnitDimension)


class PrintSettings(CreateObjectSettings):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.mmToMicronScaleFactor = 1000
        self.cellsPerUnitMM = (self.cellsPerHeatSourceDimension/self.heatSourceSizeX)*self.mmToMicronScaleFactor
        self.cellCountX = int(self.length*self.cellsPerUnitMM)
        self.cellCountY = int(self.width*self.cellsPerUnitMM)
        self.cellCountZ = int(self.thickness*self.cellsPerUnitMM)








