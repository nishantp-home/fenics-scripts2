class PrinterSpecifications:

    def __init__(self, scanSpeed=1000, laserPower=300, laserEfficiency=0.7, hatchDistance=0.05) -> None:
        self.scanSpeed = scanSpeed
        self.laserPower = laserPower
        self.laserEfficiency = laserEfficiency
        self.hatchDistance = hatchDistance