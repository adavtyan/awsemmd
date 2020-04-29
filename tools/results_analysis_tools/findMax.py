import sys

filename = sys.argv[1]
st = int(sys.argv[2])

file = open(filename, 'r')
a = file.read()
a = a.split()
max = 0
imax = -1
i=0
for ia in a[st:]:
	ia = float(ia)
	if ia>max:
		max = ia
		imax = i
	i = i + 1
print max, imax+st
file.close()
