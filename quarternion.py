import numpy

class IntegerQuarternion:
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
		return IntegerQuarternion(z, i, j, k)

	def __mul__(self, other):
		z = self.z * other.z - self.i * other.i - self.j * other.j - self.k * other.k
		i = self.i * other.z + self.z * other.i + self.j * other.k - self.k * other.j
		j = self.z * other.j + self.j * other.z - self.i * other.k + self.k * other.i
		k = self.z * other.k + self.k * other.z + self.i * other.j - self.j * other.k 
		return IntegerQuarternion(z, i, j, k)

	def conj(self):
		return IntegerQuarternion(self.z, -self.i, -self.j, -self.k)

