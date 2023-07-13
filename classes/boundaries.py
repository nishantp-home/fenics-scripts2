from fenics import *

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


class Back(SubDomain):
	tolerance = 1E-14

	def inside(self, x, on_boundary):
		return abs(x[1]) < self.tolerance and on_boundary

class Front(SubDomain):
	tolerance = 1E-14
	
	def __init__(self, specimenDimensions) -> None:
		super().__init__()
		self.specimenDimensions = specimenDimensions

	def inside(self, x, on_boundary):
		return abs(x[1] - self.specimenDimensions.width) < self.tolerance and on_boundary
	

class Top(SubDomain):
	tolerance = 1E-14

	def __init__(self, specimenDimensions) -> None:
		super().__init__()
		self.specimenDimensions = specimenDimensions

	def inside(self, x, on_boundary):
		return abs(x[2]-self.specimenDimensions.thickness) < self.tolerance and on_boundary

class Bottom(SubDomain):
	tolerance = 1E-14

	def inside(self, x, on_boundary):
		return abs(x[2]) < self.tolerance and on_boundary
	

class MYBC:
	def __init__(self, specimenDimensions, mesh, left=None, right=None, bottom=None, top=None, back=None, front=None) -> None:
		self.mesh = mesh
		self.boundary_markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim()-1).set_all(7)
		if left == "D":
			self.left = Left()			
		if right == "D":
			self.right = Right(specimenDimensions)
		if bottom == "D":
			self.bottom = Bottom()
		if top == "D":
			self.top = Top(specimenDimensions)
		if back == "D":
			self.back = Back()
		if front == "D":
			self.front = Front(specimenDimensions)
		
		self.setMarkers()


	def setMarkers(self):
		print("Assigning markers to all facets")
		self.left.mark(self.boundary_markers, 0)
		self.right.mark(self.boundary_markers, 1)
		self.bottom.mark(self.boundary_markers, 2)
		self.top.mark(self.boundary_markers, 3)
		self.back.mark(self.boundary_markers, 4)
		self.front.mark(self.boundary_markers, 5)

				