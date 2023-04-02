from fenics import *
from classes.createObjectSettings import *

class Top(SubDomain):
	tolerance = 1E-14

	def __init__(self, specimenDimensions) -> None:
		super().__init__()
		self.specimenDimensions = specimenDimensions

	def inside(self, x, on_boundary):
		return abs(x[1]-self.specimenDimensions.thickness) < self.tolerance and on_boundary

class Bottom(SubDomain):
	tolerance = 1E-14

	def inside(self, x, on_boundary):
		return abs(x[1]) < self.tolerance and on_boundary

class Left(SubDomain):
	tolerance = 1E-14

	def inside(self, x, on_boundary):
		return abs(x[0]) < self.tolerance and on_boundary

class Right(SubDomain):
	tolerance = 1E-14
	
	def __init__(self, specimenDimensions) -> None:
		super().__init__()
		self.specimenDimensions = specimenDimensions

	def inside(self, x, on_boundary):
		return abs(x[0] - self.specimenDimensions.length) < self.tolerance and on_boundary
