from abc import ABCMeta, abstractmethod

class ParameterField(metaclass = ABCMeta):
    """Interface for parameter fields"""
    
    @abstractmethod
    def __init__(self, spaceDependency=False, temperatureDependency=False, 
                 spatialExpression=None, temperatureDependentExpression=None, temperatureDependentTable=None, store=False) -> None:
           
        self.spaceDependency = spaceDependency
        self.temperatureDependency = temperatureDependency
        self.store = store
        self.spatialExpression = None
        self.temperatureDependentExpression = None
        self.temperatureDependentTable = None

        if(self.spaceDependency == True and self.temperatureDependency == True):
             raise RuntimeError("Both Space and Temperature dependent expressions for parameter fields are not supported yet.")

        if(self.spaceDependency == True and self.temperatureDependency == False):
            if(spatialExpression == None):
                raise RuntimeError("Please provide the string for spatial expression for", str(self))
            
            self.spatialExpression = spatialExpression
        
        if(self.temperatureDependency == True and self.spaceDependency == False):
            
            if(temperatureDependentExpression != None and temperatureDependentTable == None):
                self.temperatureDependentExpression = temperatureDependentExpression

            elif(temperatureDependentTable != None and temperatureDependentExpression == None):
                self.temperatureDependentTable = temperatureDependentTable
            
            else:
                raise RuntimeError("Please provide the string for temperature dependent expression or csv file containing temperature dependent table for", str(self))

        
