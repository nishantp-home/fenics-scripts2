class SimulationParameters:

    __timeToleranceFactor = 1e-6

    def __init__(self, storedFrameCount=5, totalTime=0.001, trapezoidalParameter=1.0) -> None:
        self.trapezoidalParameter = trapezoidalParameter
        self.storedFrameCount = storedFrameCount
        self.totalTime = totalTime
        self.dt = 1E-5
        self.timeTolerance = self.__timeToleranceFactor*self.dt
        self.timeStepCount = int(self.totalTime/self.dt)
        self.storeFrequency = int(self.timeStepCount/self.storedFrameCount)
        

class PrintSimulationParameters:

    __timeToleranceFactor = 1e-6

    def __init__(self, printSettings, storedFrameCount=5, trapezoidalParameter=1.0) -> None:
        self.printSettings = printSettings
        self.trapezoidalParameter = trapezoidalParameter
        self.storedFrameCount = storedFrameCount
        self.totalTime = self.printSettings.length/self.printSettings.scanSpeed        
        self.dt = 1E-5
        self.timeTolerance = self.__timeToleranceFactor*self.dt
        self.timeStepCount = int(self.totalTime/self.dt)
        self.storeFrequency = int(self.timeStepCount/self.storedFrameCount)