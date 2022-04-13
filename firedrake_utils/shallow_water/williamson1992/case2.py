
from math import pi

import firedrake as fd
from firedrake_utils.planets import earth
from firedrake_utils import units

### === --- constants --- === ###

# gravitational constant * reference depth
gh0 = 2.94e4*pow(units.metre/units.second,2)

# reference depth
h0 = gh0/earth.gravity

# days taken for velocity to travel circumference
period = 12.0

# reference velocity
u0 = 2*pi*earth.radius/(period*earth.day)

### === --- analytical solution --- === ###

# coriolis parameter f
def coriolis_expression( x,y,z ):
    return 2*earth.omega*z/earth.radius

def coriolis_function( x,y,z, Vf, name="coriolis" ):
    return fd.Function( Vf, name=name ).interpolate(coriolis_expression(x,y,z))

# velocity field u
def velocity_expression( x,y,z ):
    return fd.as_vector([ -u0*y/earth.radius, u0*x/earth.radius, 0.0 ])

def velocity_function( x,y,z, V1, name="velocity" ):
    return fd.Function( V1, name=name ).project(velocity_expression(x,y,z))

# depth field h
def depth_expression( x,y,z ):
    return h0 - ( earth.radius*earth.omega*u0 + 0.5*u0*u0 )*(z*z/(earth.radius*earth.radius))/earth.gravity

def depth_function( x,y,z, V2, name="depth" ):
    return fd.Function( V2, name=name ).interpolate(depth_expression(x,y,z))

