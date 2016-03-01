from quarternion import IntegerQuarternion
from finitefield.finitefield import FiniteField

def main():
	x = IntegerQuarternion(1, 2, 3, 4)
	print x
	print x + x
	print x * x
	F5 = FiniteField(5,1)
	x = F5(2)
	print x * x

if __name__ == '__main__':
	main()