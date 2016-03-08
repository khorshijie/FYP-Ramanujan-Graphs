import numpy
from finitefield.finitefield import FiniteField
from finitefield.polynomial import *
from finitefield.modp import *
from igraph import *

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
		k = self.z * other.k + self.k * other.z + self.i * other.j - self.j * other.i
		return Quarternion(z, i, j, k)

	def __eq__(self, other):
		return isinstance(other, Quarternion) and self.z == other.z and self.i == other.i and self.j == other.j and self.k == other.k 

	def conj(self):
		return Quarternion(self.z, -self.i, -self.j, -self.k)

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

# TODO: Improve code quality here
def get_distinguished_set(p):
	Sp = []
	for i in range(0, p):
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

def get_spq(p, q):
	Sp = get_distinguished_set(p)
	Fq = FiniteField(q,1)
	FFQq = FiniteFieldQuarternion(Fq)
	Spq = []
	for quar in Sp:
		Spq.append(embed_into_pgl2p(FFQq.convert_to_ffq(quar).to_matrix(), q))
	return Spq

def get_gl2p(p):
	gl2p = []
	FFp = FiniteField(p)
	for i in range(p):
		for j in range(p):
			for k in range(p):
				for l in range(p):
					if (i * l - k * j) % p != 0:
						gl2p.append(numpy.matrix([[FFp(i), FFp(j)], [FFp(k),FFp(l)]]))
	return gl2p

def get_sl2p(p):
	gl2p = []
	FFp = FiniteField(p)
	for i in range(p):
		for j in range(p):
			for k in range(p):
				for l in range(p):
					if (i * l - k * j) % p == 1:
						gl2p.append(numpy.matrix([[FFp(i), FFp(j)], [FFp(k),FFp(l)]]))
	return gl2p

def get_psl2p(p):
	sl2p = get_sl2p(p)
	return [x for x in sl2p if numpy.array_equal(x, embed_into_psl2p(x, p))]

def get_pgl2p(p):
	gl2p = get_gl2p(p)
	return [x for x in gl2p if (numpy.array_equal(x,embed_into_pgl2p(x, p)))]

def embed_into_pgl2p(mat, p):
	FFp = FiniteField(p)
	inv = FFp(1)
	if mat[0,0] != FFp(0):
		inv = inv / mat[0,0]		
	else:
		inv = inv / mat[0,1]
	diag = numpy.matrix([[inv, FFp(0)], [FFp(0), inv]])
	return diag * mat	

def embed_into_psl2p(mat, p):
	FFp = FiniteField(p)
	diag = numpy.matrix([[FFp(-1), FFp(0)], [FFp(0), FFp(-1)]])
	if mat[0,0] != FFp(0):
		if mat[0,0] > FFp((p-1)/2):
			return diag * mat
		else :
			return mat
	else : 
		if mat[0,1] > FFp((p-1)/2):
			return diag * mat
		else :
			return mat

def legendre(p,q):
	if p % q == 0:
		return 0
	elif p == 1:
		return 1
	elif p > q:
		return legendre(p % q, q)
	else:
		exp = (p-1)*(q-1)/4
		if exp %2 == 1:
			return legendre(q, p)
		else : 
			return -legendre(q,p)

def draw_Xpq(p, q):
	Spq = get_spq(p,q)
	mat_to_int = {}
	int_to_mat = {}
	if(legendre(p, q) == 1):
		psl2q = get_psl2p(q)
		edges = []
		for x in xrange(len(psl2q)):
			for op in Spq:
				product = embed_into_psl2p(psl2q[x] * op, q)
				edges.append((x, find_matrix_list(psl2q, product)))
		g = Graph()
		g.add_vertices(len(psl2q))
		g.add_edges(edges)
		layout = g.layout("drl")
		plot(g, layout = layout)
	else :
		pgl2q = get_pgl2p(q)
		edges = []
		for x in xrange(len(pgl2q)):
			for op in Spq:
				product = embed_into_pgl2p(pgl2q[x] * op, q)
				edges.append((x, find_matrix_list(pgl2q, product)))
		g = Graph()
		g.add_vertices(len(pgl2q))
		g.add_edges(edges)
		layout = g.layout("drl")
		plot(g, layout = layout)

def find_matrix_list(lst, x):
	cnt = 0
	for item in lst:
		if numpy.array_equal(item, x):
			return cnt
		else:
			cnt = cnt + 1
	return cnt
