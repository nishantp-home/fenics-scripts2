from .parameterFields import ParameterField


class HeatConductivity(ParameterField):
    def __init__(self, constant=1.0, spaceDependency=False, temperatureDependency=False, 
                 spatialExpression=None, temperatureDependentExpression=None, temperatureDependentTable=None, store=False) -> None:
        self.constant = constant
        super().__init__(spaceDependency, temperatureDependency, spatialExpression, temperatureDependentExpression, temperatureDependentTable, store)
        