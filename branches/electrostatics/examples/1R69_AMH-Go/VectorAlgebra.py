from math import *

def det3(a):
    res = a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
    res += a[0][1]*(a[1][2]*a[2][0]-a[2][2]*a[1][0])
    res += a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])
    return res

def sgn(a):
    if a>0: return 1
    elif a<0: return -1
    else: return 0

def vector(p1, p2):
    return [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]

def vproduct(a, b):
    if type(a)==type([]) and type(b)==type([]):
        return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
    elif type(b)==type([]):
        return [a*b[0], a*b[1], a*b[2]]
    elif type(a)==type([]):
        return [a[0]*b, a[1]*b, a[2]*b]
    return a*b

def vcross_product(a, b):
    cx = a[1]*b[2]-a[2]*b[1]
    cy = a[2]*b[0]-a[0]*b[2]
    cz = a[0]*b[1]-a[1]*b[0]
    return [cx, cy, cz];

def vabs(a):
    return sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2))

def vangle(a, b):
    return acos(vproduct(a, b)/(vabs(a)*vabs(b)))

def dihedral_angle(v1, v2, v3):
    n1 = vcross_product(v1, v2)
    n2 = vcross_product(v2, v3)
    y = vproduct( vproduct(vabs(v2), v1), n2 )
    x = vproduct( n1, n2 )
    return 180*atan2(y, x)/pi

def mproduct(a, b):
    na = len(a)
    if na>0: ma = len(a[0])
    else: return []
    for ia in a:
        if len(ia)!=ma: return []
    
    nb = len(b)
    if nb>0: mb = len(b[0])
    else: return []
    for ib in b:
        if len(ib)!=mb: return []

    if ma!=nb: return []

    c=[]
    for i in range(0,na):
        c.append([])
        for j in range(0,mb):
            cij = 0
            for k in range(0,ma):
                cij += a[i][k]*b[k][j]
            c[i].append(cij)

    return c


# phi rotate XYZ around Z until X coincides with K
# theta rotate around K until Z coincideswith Z'
# psi rotate around Z' until K coincides with X'
def R(phi, theta, psi):
    R11 = cos(phi)*cos(psi)-sin(phi)*cos(theta)*sin(psi)
    R12 = sin(phi)*cos(psi)+cos(phi)*cos(theta)*sin(psi)
    R13 = sin(theta)*sin(psi)

    R21 = -cos(phi)*sin(psi)-sin(phi)*cos(theta)*cos(psi)
    R22 = -sin(phi)*sin(psi)+cos(phi)*cos(theta)*cos(psi)
    R23 = sin(theta)*cos(psi)

    R31 = sin(phi)*sin(theta)
    R32 = -cos(phi)*sin(theta)
    R33 = cos(theta)
    
    return [[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]]

