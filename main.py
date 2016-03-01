from quarternion import *
# Finite field library obtained from Jeremy Kun's blog
# https://github.com/j2kun/finite-fields
from finitefield.finitefield import *
import numpy

def main():
	# testing quarternion
	x = Quarternion(1, 2, 3, 4)
	print x
	print x + x
	print x * x
	
	# testing finite field quarternion
	F5 = FiniteField(5, 1)
	FFQ5 = FiniteFieldQuarternion(F5)
	y = FFQ5(F5(1), F5(2), F5(3), F5(4))
	print y.to_matrix()
	sp = get_distinguished_set(5)
	for i in sp:
		print i

if __name__ == '__main__':
	main()