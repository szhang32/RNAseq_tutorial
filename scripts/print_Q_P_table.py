
print("Char\tInt\tQ\tP")
for Q in range(0, 43):
	#if (Q == 14 or Q == 28):
	#	print("Char\tInt\tQ\tP")
	ascii = Q + 33
	print(chr(ascii), end = '\t')	
	print(ascii, end = '\t')
	print(Q, end = '\t')
	P = 10 ** (-1 * Q/10)
	print('%.6f' % P, end = '\t')
	print()
