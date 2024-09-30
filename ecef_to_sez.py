# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
# Finds ECEF vector to the origin of the SEZ reference frame
# Parameters:
# o_x_km, o_y_km, o_z_km : ECEF origin of SEZ frame
# x_km, y_km, z_km : ECEF position
# Output:
# s_km, e_km, z_km
#
# Written by Anna Kosnic
#
# import Python modules
import sys # argv
import math as m

# "constants"
R_E_KM = 6378.1363
E_E = 0.081819221456

# helper functions
def calc_denom(ecc, lat):
    return m.sqrt(1-ecc**2 *(m.sin(lat))**2)

# initialize script arguments
o_x_km  = float('nan')
o_y_km  = float('nan')
o_z_km  = float('nan')
x_km    = float('nan')
y_km    = float('nan')
z_km    = float('nan')

# parse script arguments
if len(sys.argv)==7:
    o_x_km = float(sys.argv[1])
    o_y_km = float(sys.argv[2])
    o_z_km = float(sys.argv[3])
    x_km   = float(sys.argv[4])
    y_km   = float(sys.argv[5])
    z_km   = float(sys.argv[6])
else:
    print(\
        'Usage: '\
        'python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km'\
            )
    exit()


# write script below this line
P_ECEF = (x_km, y_km, z_km) # ECEF position vector 
O_ECEF   = (o_x_km, o_y_km, o_z_km) # ECEF origin of SEZ vector

T = [P_ECEF[0]-O_ECEF[0], P_ECEF[1]-O_ECEF[1], P_ECEF[2]-O_ECEF[2]] # translation


# ECEF to LLH of ground site (SEZ origin)
# calculate longitude
lon_rad = m.atan2(o_y_km,o_x_km)
lon_deg = lon_rad * 180/m.pi

# inital guess
lat_rad = m.asin(o_z_km/m.sqrt(o_x_km**2 + o_y_km**2 + o_z_km**2))
r_lon_km = m.sqrt(o_x_km**2 + o_y_km**2)
prev_lat_rad = float('nan')

# iterative solver
count = 0
c_E = float('nan')

while (m.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad>10e-12)) and count<10:
    denom = calc_denom(E_E, lat_rad)
    c_E = R_E_KM/denom

    prev_lat_rad = lat_rad
    lat_rad = m.atan((o_z_km + c_E*E_E**2 * m.sin(lat_rad)) / r_lon_km)

    count = count + 1

hae_km = r_lon_km/m.cos(lat_rad)-c_E

# applying inverse rotations
Ry = [m.cos(lon_rad)*T[0] + m.sin(lon_rad)*T[1], -m.sin(lon_rad)*T[0] + m.cos(lon_rad)*T[1], T[2]]
Rz = [m.sin(lat_rad)*Ry[0] + -m.cos(lat_rad)*Ry[2], Ry[1], m.cos(lat_rad)*Ry[0] + m.sin(lat_rad)*Ry[2]]

s_km = Rz[0]
e_km = Rz[1]
z_km = Rz[2]

print(s_km)
print(e_km)
print(z_km)