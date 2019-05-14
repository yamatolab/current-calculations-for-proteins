from __future__ import print_function
import os, sys
import numpy

topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path:
    sys.path.insert(0, topdir)
import clog as logger

def cal_com_crd_vel(crd_or_vel, masses):
    """Calculate the coordinate and velocity of the center of mass."""
    
    # com : center of mass
    mass_com = numpy.sum(masses)
    com_crd_vel = numpy.sum(masses[:, None] * crd_or_vel, 0) / mass_com
    return com_crd_vel

    # above code means below code.
    # x_com = numpy.sum(masses * crd_or_vel[:, 0])
    # y_com = numpy.sum(masses * crd_or_vel[:, 1])
    # z_com = numpy.sum(masses * crd_or_vel[:, 2])
    # return numpy.array([x_com, y_com, z_com]) / mass_com

def cal_rotation_omega(crd, vel, masses):
    """Calculate the velocity of rotation coordinate."""

    # calculate the total angular momentum : L[3]
    r_x_mv = masses[:, None] * numpy.cross(crd, vel)
    L = numpy.sum( r_x_mv, 0 )

    # calculate the inertia tensor : I[3,3]
    xs, ys, zs = crd[:, 0], crd[:, 1], crd[:, 2]

    I_xx = numpy.sum(masses * (ys*ys + zs*zs), 0)
    I_yy = numpy.sum(masses * (zs*zs + xs*xs), 0)
    I_zz = numpy.sum(masses * (xs*xs + ys*ys), 0)
    I_xy = - numpy.sum(masses * xs * ys, 0)
    I_xz = - numpy.sum(masses * xs * zs, 0)
    I_yx = I_xy
    I_yz = - numpy.sum(masses * ys * zs, 0)
    I_zx = I_xz
    I_zy = I_yz

    I = numpy.array([
        [I_xx, I_xy, I_xz],
        [I_yx, I_yy, I_yz],
        [I_zx, I_zy, I_zz]
        ])

    # angular momentum: omega = I^-1 L
    I_inv = numpy.linalg.inv(I)
    omegas = numpy.dot(I_inv, L)
    return omegas

def get_crd_vel_trans_rot_removed(masses, crd, vel, target_atoms,
                                rem_trans, rem_rotate):
    target_indices = numpy.array(target_atoms) - 1
    target_masses = numpy.array(masses)[target_indices]

    # get the crd and vel in target atoms
    target_crd = crd[target_indices]
    target_vel = vel[target_indices]

    logger.debug('*** Translation and Rotation ***')
    logger.debug('    calculating within target atoms...\n')
    msg = '{:>30} = ({:8.5f}, {:8.5f}, {:8.5f}) [{:}]'
    # msg = '{:>30} = ({:10.7e}, {:10.7e}, {:10.7e}) [{:}]'

    if rem_trans:
        # get the center of mass for the coordinate
        target_r_com = cal_com_crd_vel(target_crd, target_masses)
        x_com, y_com, z_com = [float(x) for x in target_r_com] # ndarray => float
        logger.info_cycle(msg.format('crd of center of mass', 
                x_com, y_com, z_com, 'A'))

        # get the center of mass for the velocity
        target_v_com = cal_com_crd_vel(target_vel, target_masses)
        vx_com, vy_com, vz_com = [float(v) for v in 1000.0*target_v_com] # ndarray => float
        logger.info_cycle(msg.format('vel of center of mass', 
                vx_com, vy_com, vz_com, 'A/ps'))

        new_crd = crd - target_r_com
        new_vel = vel - target_v_com
        target_crd_minus_com = target_crd - target_r_com
        target_vel_minus_com = target_vel - target_v_com

    else:
        logger.info_cycle(20*' ' + 'remove_trans is off.')
        new_crd = crd
        new_vel = vel
        target_crd_minus_com = target_crd
        target_vel_minus_com = target_vel

    if rem_rotate:

        target_omegas = cal_rotation_omega(
                target_crd_minus_com, target_vel_minus_com, target_masses)

        # log the omega values.
        x_omega, y_omega, z_omega = [float(v) for v in 1000.0*target_omegas]
        logger.info_cycle(msg.format('angular vel for rotation', 
                x_omega, y_omega, z_omega, 'rad/ps'))

        # calculate rotational motion for all of the system using target_omegas.
        target_v_rot = numpy.cross(target_omegas, new_crd)
        new_vel = new_vel - target_v_rot

    else:
        logger.info_cycle(20*' ' + 'remove_rotate is off.')

    logger.info_cycle()

    return new_crd, new_vel

if __name__ == '__main__':

    crd = numpy.array([[-1.2,0.5,-1.3],[0.8,-1.2,1.1]])
    vel = numpy.array([[ 1.0,0.3,0.9],[-0.7,0.9,1.5]])
    mass = numpy.array([1.1, 0.9])

    # cal_rotation_omega(crd, vel, mass)

    new_crd, new_vel = get_crd_vel_trans_rot_removed(
            mass, crd, vel, [1,2], True, True)

    print(new_vel)
