
from math import pi

import firedrake as fd
from firedrake_utils.planets import earth
from firedrake_utils import units

from firedrake_utils.shallow_water.williamson1992 import case2

# # # === --- constants --- === # # #

# reference depth
h0 = 5960*units.metre

# reference velocity
u0 = 20*units.metre/units.second

# mountain parameters
mountain_height = 2000*units.metre

mountain_radius = pi/9.

# different lambda_c because atan_2 used for angle
mountain_centre_lambda = -pi/2.
mountain_centre_theta = pi/6.

# # # === --- analytical solution --- === # # #


# coriolis parameter f
def coriolis_expression(x, y, z):
    return case2.coriolis_expression(x, y, z)


def coriolis_function(x, y, z, Vf, name="coriolis"):
    return fd.Function(Vf, name=name).interpolate(coriolis_expression(x, y, z))


# velocity field u
def velocity_expression(x, y, z, uref=u0):
    return case2.velocity_expression(x, y, z, uref=uref)


def velocity_function(x, y, z, V1, uref=u0, name="velocity"):
    v = fd.Function(V1, name=name)
    return v.project(velocity_expression(x, y, z, uref=uref))


# depth field h
def depth_expression(x, y, z,
                     radius=mountain_radius,
                     height=mountain_height,
                     theta_c=mountain_centre_theta,
                     lambda_c=mountain_centre_lambda):

    lambda_x = fd.atan_2(y/earth.radius, x/earth.radius)
    theta_x = fd.asin(z/earth.radius)

    radius2 = pow(radius, 2)
    lambda2 = pow(lambda_x - lambda_c, 2)
    theta2 = pow(theta_x - theta_c,  2)

    min_arg = fd.Min(radius2, theta2 + lambda2)

    return height*(1 - fd.sqrt(min_arg)/radius)


def depth_function(x, y, z, V2,
                   radius=mountain_radius,
                   height=mountain_height,
                   theta_c=mountain_centre_theta,
                   lambda_c=mountain_centre_lambda,
                   name="depth"):

    h = fd.Function(V2, name=name)
    hexp = depth_expression(x, y, z,
                            radius=radius,
                            height=height,
                            theta_c=theta_c,
                            lambda_c=lambda_c)
    return h.interpolate(hexp)
