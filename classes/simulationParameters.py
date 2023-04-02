class SimulationParameters:

    _timeToleranceFactor = 1e-6

    def __init__(self, storedFrameCount=5, totalTime=0.1, trapezoidalParameter=1.0) -> None:
        self.trapezoidalParameter = trapezoidalParameter
        self.storedFrameCount = storedFrameCount
        self.totalTime = totalTime
        self.dt = 1E-5
        self.timeTolerance = self._timeToleranceFactor*self.dt
        self.timeStepCount = int(self.totalTime/self.dt)
        self.storeFrequency = int(self.timeStepCount/self.storedFrameCount)
        