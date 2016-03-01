import numpy
from finitefield.finitefield import FiniteField
from finitefield.polynomial import *
from finitefield.modp import *

class Quarternion():
	def __init__(self, z, i, j, k):
		self.z = z
		self.i = i
		self.j = j
		self.k = k

	def __str__(self):
		return str(self.z) + " + " + str(self.i) + "i + " + str(self.j) + "j + " + str(self.k) + "k"

	def __add__(self, other):
		z = self.z + other.z
		i = self.i + other.i
		j = self.j + other.j
		k = self.k + other.k
		return Quarternion(z, i, j, k)

	def __mul__(self, other):
		z = self.z * other.z - self.i * other.i - self.j * other.j - self.k * other.k
		i = self.i * other.z + self.z * other.i + self.j * other.k - self.k * other.j
		j = self.z * other.j + self.j * other.z - self.i * other.k + self.k * other.i
		k = self.z * other.k + self.k * other.z + self.i * other.j - self.j * other.k 
		return Quarternion(z, i, j, k)

	def __eq__(self, other):
		return isinstance(other, Quarternion) and self.z == other.z and self.i == other.i and self.j == other.j and self.k == other.k 

	def conj(self):
		return Quarternion(self.z, -self.i, -self.j, -self.k)

# TODO: Improve code quality here
def get_distinguished_set(p):
	Sp = []
	for i in range(-p, p):
		for j in range(-p, p):
			for k in range(-p, p):
				for l in range(-p, p):
					if i ** 2 + j ** 2 + k ** 2 + l ** 2 == p:
						if i % 2 == 1 and j % 2 == 0 and k % 2 == 0 and l % 2 == 0 and i > 0:
							Sp.append(Quarternion(i, j, k, l))
						elif i % 2 == 0 and j % 2 == 1 and k % 2 == 1 and l % 2 == 1:
							if i != 0:
								Sp.append(Quarternion(i, j, k, l))
							elif (j != 0 and j > 0) or (k != 0 and k > 0) or (l != 0 and l > 0):
								Sp.append(Quarternion(i, j, k, l)) 
	return Sp

def FiniteFieldQuarternion(FiniteField):

	if hasattr(FiniteField, 'primeSubfield'):
		primeSubfield = FiniteField.primeSubfield
	else:
		primeSubfield = FiniteField

	p = primeSubfield.p

	def solve_for_xy():		
		for i in range(p):
			for j in range(p):
				if i * i + j * j + primeSubfield(1) == primeSubfield(0):
					return i, j
		raise Exception("Your field is bad.")

	x, y = solve_for_xy()

	class FFQ(Quarternion):
		@classmethod
		def convert_to_ffq(self,q):
			return FFQ(FiniteField(q.z), FiniteField(q.i), FiniteField(q.j), FiniteField(q.k))

		def to_matrix(self):
			return numpy.matrix([[self.z + self.i * x + self.k * y, -self.i * y + self.j + self.k * x], 
				[-self.i * y - self.j + self.k * x, self.z - self.i * x - self.k * y]], dtype = FiniteField)

	return FFQ