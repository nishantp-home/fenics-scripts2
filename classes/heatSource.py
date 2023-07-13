from abc import ABCMeta, abstractmethod

class IHeatSource(metaclass = ABCMeta):
    """Interface for parameter fields"""

    @abstractmethod
    def __init__(self) -> None:
        pass


class CubicHeatSource(IHeatSource):
    """Cubic heat Source specification"""
    def __init__(self, cubeEdgeLength = 0.1) -> None:
        super().__init__()
        self.cubeEdgeLength = cubeEdgeLength

        