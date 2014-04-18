import sys

filename = sys.argv[1]
st = int(sys.argv[2])

file = open(filename, 'r')
a = file.read()
a = a.split()
min = 0
imin = -1
i=0
for ia in a[st:]:
	ia = float(ia)
	if ia<min or imin==-1:
		min = ia
		imin = i
	i = i + 1
print min, imin+st
file.close()
