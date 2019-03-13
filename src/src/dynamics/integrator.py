from __future__ import print_function

# standard modules
import os, sys
import time
import numpy

# curp modules
topdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if topdir not in sys.path: sys.path.insert(0, topdir)
from utility import TimeStore
import clog as logger

class NEVIntegrator(TimeStore):

    coef = 4.184*10.0**(-4) # dt*f/m (fs*kcal/mol/A/u) => v(A/fs)
    gas_const = 8.3144621 # gas constant (J/K/mol)

    def __init__(self, topology, setting, interact_table, crd=None, vel=None):
        TimeStore.__init__(self)
        self.setup(topology, setting, interact_table)
        self.__crd = crd
        self.__vel = vel

    def setup(self, topology, setting, interact_table):
        self.__tpl = topology
        self.__setting = setting

        # decide force calculator
        import twobody
        TwoBodyCalculator = twobody.get_calculator(setting.curp.potential)
        self.__tbf = TwoBodyCalculator(topology, setting)
        self.__tbf.setup(interact_table, check=False)

        self.__interact_table = [ numpy.array(t) for t in interact_table ] 

        # decide integrator
        params = setting.dynamics
        if params.integrator == 'vverlet':
            self.integrate = self.do_vverlet
        elif params.integrator == 'leapfrog':
            self.integrate = self.do_leapfrog
        else:
            pass

    def run(self, data=None):
        """Run the one step integrator."""

        if data:
            cstep, (crd, vel, pbc) = data
        else:
            crd, vel = self.__crd, self.__vel

        logger.debug('    calculate dynamics ...')

        masses   = self.__tpl.get_atom_info()['masses']
        force = self.cal_force(crd)
        params = self.__setting.dynamics

        results = self.integrate(crd, force, vel, masses, params)

        return results

    def cal_force(self, crd):
        return self.__tbf.cal_force(crd)

    def integrate(self, crd, frc, vel, masses):
        pass
    
    def do_vverlet(self, crd, frc, vel, masses, params):
        """Integrate the coordinate and velocity according to
        velocity verlet algorithm.
        """

        # istp     : current step
        # crd      : current step's coordinate
        # vel      : current step's velocity
        # frc      : current step's force
        # crd_next : next step's coordinate
        # vel_next : next step's velocity
        # frc_next : next step's force

        # Preparation of the dynamcs on 0th step
        ms    = numpy.array(masses)[:, None]
        dt    = params.dt * 1000.0
        nstep = params.num_steps
        coef  = self.coef

        # log the information of the restart file
        logger.info('*** The information of the restart file ***')
        # log temperature
        ek, temp = self.cal_temp(vel, ms)
        self.output_temp(ek, temp)
        # log energy
        if self.__setting.output.output_energy: self.__tbf.output_energy()
        # log forces
        if logger.is_debug(): self.__tbf.output_force()

        # do dynamics
        for istp_1 in range(nstep):
            istp = istp_1+1

            t0 = time.time()

            logger.set_curstep(istp)
            logger.info_cycle('*** ISTEP = {}, CURRENT TIME = {} ***'
                .format(istp, istp * dt))

            # calculate next step's coordinate
            # r(t+dt) = r(t) + dt*v(t) + dt^2/2m * F(t)
            # crd_next = crd + dt*vel + 0.5*dt*dt * frc / ms * coef
            crd_next = crd + dt*vel + 0.5*dt*dt * frc / ms * coef
            
            # calculate next step's forces from coordinate
            frc_next = self.cal_force(crd_next)

            # calculate next step's velocity
            # v(t+dt) = v(t) + dt/2m * { F(t+dt) + F(t) }
            vel_next = vel + 0.5*dt* (frc_next+frc) / ms * coef

            self.store_time('integrator', time.time()-t0)
            yield istp, (crd_next, vel_next)

            crd = crd_next
            frc = frc_next
            vel = vel_next

            # log temperature
            ek, temp = self.cal_temp(vel, ms)
            self.output_temp(ek, temp)
            # log energy
            if self.__setting.output.output_energy: self.__tbf.output_energy()
            # log forces
            if logger.is_debug(): self.__tbf.output_force()

    def do_leapfrog(self, crd, frc, vel, masses, params):
        """Integrate the coordinate and velocity according to
        leap frog algorithm.
        """

        # istp     : current step(t)
        # crd      : current step(t)'s coordinate
        # frc      : current step(t)'s force
        # vel      : current step(t-dt/2)'s velocity
        # crd_next : next step(t+dt)'s coordinate
        # vel_next : next step(t+dt/2))'s velocity
        # frc_next : next step(t+dt)'s force
        # vel_dt = 0.5*(vel + vel_next)

        # Preparation of the dynamcs on 0th step
        ms    = numpy.array(masses)[:, None]
        # dt    = params.dt
        # dt    = params.dts
        dt = params.dt * 1000.0
        nstep = params.num_steps
        coef  = self.coef

        # log the information of the restart file
        logger.info('*** The information of the restart file ***')
        logger.info("    Velocity is 0+1/2-th step' informations")
        # log temperature
        ek, temp = self.cal_temp(vel, ms)
        self.output_temp(ek, temp)
        # log energy
        if self.__setting.output.output_energy: self.__tbf.output_energy()
        # log forces
        if logger.is_debug(): self.__tbf.output_force()

        # do dynamics
        for istp_1 in range(nstep):
            istp = istp_1+1

            t0 = time.time()

            logger.set_curstep(istp)
            logger.info_cycle('*** ISTEP = {}, CURRENT TIME = {} ***'
                .format(istp, istp * dt))

            # calculate next step's velocity
            # v(t+dt/2) = v(t-dt/2) + dt * F(t) /m
            vel_next = vel + dt*frc/ms * coef
            vel_t   = 0.5*( vel + vel_next )

            # calculate next step's coordinate
            # r(t+dt) = r(t) + dt*v(t+dt/2)
            crd_next = crd + dt*vel_next
            
            # calculate next step's forces from coordinate
            frc_next = self.cal_force(crd_next)

            self.store_time('integrator', time.time()-t0)
            yield istp, (crd_next, vel_next)

            crd = crd_next
            frc = frc_next
            vel = vel_next

            # log temperature
            ek, temp = self.cal_temp(vel_t, ms)
            self.output_temp(ek, temp)
            # log energy
            if self.__setting.output.output_energy: self.__tbf.output_energy()
            # log forces
            if logger.is_debug(): self.__tbf.output_force()

    def cal_temp(self, vel, ms):
        """Calculate the temperature and Kenitic energy."""
        gas  = self.gas_const
        dof  = vel.size # defree of freedom

        energy_tmp = numpy.sum( ms * vel**2)

        energy = 0.5 * energy_tmp / self.coef
        temp   = 10.0**7 * energy_tmp / dof / gas
        return energy, temp

    def output_temp(self, energy, temp):
        """Output the temperature and Kinetic energy into log file."""
        logger.info_cycle('    Temperature    : {:>8.2f} (K) '.format(temp))
        logger.info_cycle('    Kinetic energy : {:>8.5f} (K) '.format(energy))

