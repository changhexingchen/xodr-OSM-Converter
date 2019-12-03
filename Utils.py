import math
from geompreds import orient2d, incircle
from opendrivepy.point import Point

a = 6378137
b = 6356752.3142
f = (a - b) / a
e_sq = f * (2-f)

def point_distance(pointa, pointb):
    return math.sqrt((pointa.x - pointb.x) ** 2 + (pointa.y - pointb.y) ** 2)

def line(p1, p2):
    A = (p1.y - p2.y)
    B = (p2.x - p1.x)
    C = (p1.x*p2.y - p2.x*p1.y)
    return A, B, -C

def intersection(L1, L2):
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return Point(x,y)
    else:
        return False

def line_cross(line1, line2):
    L1 = line(line1[0], line1[1])
    L2 = line(line2[0], line2[1])

    R = intersection(L1, L2)
    if R:
        return R
    else:
        x1 = line1[0].x
        x2 = line1[1].x
        x3 = line2[0].x
        x4 = line2[1].x

        y1 = line1[0].y
        y2 = line1[1].y
        y3 = line2[0].y
        y4 = line2[1].y
        return Point((x1+x2+x3+x4)/4,(y1+y2+y3+y4)/4)

def orient_node(nodea, nodeb, nodec):
    return orient2d((nodea.x, nodea.y), (nodeb.x, nodeb.y), (nodec.x, nodec.y))


def find_diagonal(nodes):
    # find the diagonal node of node[0]
    if orient_node(nodes[0], nodes[1], nodes[2]) * orient_node(nodes[3], nodes[1], nodes[2]) < 0:
        return 3
    elif orient_node(nodes[0], nodes[1], nodes[3]) * orient_node(nodes[2], nodes[1], nodes[3]) < 0:
        return 2
    else:
        return 1




def enu_to_ecef(xEast, yNorth, zUp, lat0, lon0, h0):
    lamb = math.radians(lat0)
    phi = math.radians(lon0)
    s = math.sin(lamb)
    N = a / math.sqrt(1 - e_sq * s * s)

    sin_lambda = math.sin(lamb)
    cos_lambda = math.cos(lamb)
    sin_phi = math.sin(phi)
    cos_phi = math.cos(phi)

    x0 = (h0 + N) * cos_lambda * cos_phi
    y0 = (h0 + N) * cos_lambda * sin_phi
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda

    t = cos_lambda * zUp - sin_lambda * yNorth

    zd = sin_lambda * zUp + cos_lambda * yNorth
    xd = cos_phi * t - sin_phi * xEast 
    yd = sin_phi * t + cos_phi * xEast

    x = xd + x0 
    y = yd + y0 
    z = zd + z0 

    return x, y, z

def ecef_to_geodetic(x, y, z):
   # Convert from ECEF cartesian coordinates to 
   # latitude, longitude and height.  WGS-84
    x2 = x ** 2 
    y2 = y ** 2 
    z2 = z ** 2 

    a = 6378137.0000    # earth radius in meters
    b = 6356752.3142    # earth semiminor in meters 
    e = math.sqrt (1-(b/a)**2) 
    b2 = b*b 
    e2 = e ** 2 
    ep = e*(a/b) 
    r = math.sqrt(x2+y2) 
    r2 = r*r 
    E2 = a ** 2 - b ** 2 
    F = 54*b2*z2 
    G = r2 + (1-e2)*z2 - e2*E2 
    c = (e2*e2*F*r2)/(G*G*G) 
    s = ( 1 + c + math.sqrt(c*c + 2*c) )**(1/3) 
    P = F / (3 * (s+1/s+1)**2 * G*G) 
    Q = math.sqrt(1+2*e2*e2*P) 
    ro = -(P*e2*r)/(1+Q) + math.sqrt((a*a/2)*(1+1/Q) - (P*(1-e2)*z2)/(Q*(1+Q)) - P*r2/2) 
    tmp = (r - e2*ro) ** 2 
    U = math.sqrt( tmp + z2 ) 
    V = math.sqrt( tmp + (1-e2)*z2 ) 
    zo = (b2*z)/(a*V) 

    height = U*( 1 - b2/(a*V) ) 
    
    lat = math.atan( (z + ep*ep*zo)/r ) 

    temp = math.atan(y/x) 
    if x >=0 :    
        long = temp 
    elif (x < 0) & (y >= 0):
        long = math.pi + temp 
    else :
        long = temp - math.pi 

    lat0 = lat/(math.pi/180) 
    lon0 = long/(math.pi/180) 
    h0 = height 

    return lat0, lon0, h0

def enu_to_geodetic(xEast, yNorth, zUp, lat_ref, lon_ref, h_ref):

    x,y,z = enu_to_ecef(xEast, yNorth, zUp, lat_ref, lon_ref, h_ref)

    return ecef_to_geodetic(x,y,z)