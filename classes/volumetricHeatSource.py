from .parameterFields import ParameterField


class VolumetricHeatSource(ParameterField):
    def __init__(self, constant=0.0, spaceDependency=False, temperatureDependency=False, 
                 spatialExpression=None, temperatureDependentExpression=None, temperatureDependentTable=None, store=False) -> None:
        self.constant = constant
        super().__init__(spaceDependency, temperatureDependency, spatialExpression, temperatureDependentExpression, temperatureDependentTable, store)


class MovingVolumetricHeatSource(VolumetricHeatSource):
    def __init__(self, printSettings, constant=0, spaceDependency=False, temperatureDependency=False, spatialExpression=None, temperatureDependentExpression=None, temperatureDependentTable=None, store=False) -> None:
        super().__init__(constant, spaceDependency, temperatureDependency, spatialExpression, temperatureDependentExpression, temperatureDependentTable, store)

        self.printSettings = printSettings

