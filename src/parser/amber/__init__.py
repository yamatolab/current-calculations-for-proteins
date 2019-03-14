
import trajectory
coordinate_dict = {'ascii' : trajectory.CoordinateParser,
                   'netcdf': trajectory.NetCDFCoordinateReader}

velocity_dict   = {'ascii' : trajectory.VelocityParser,
                   'netcdf': trajectory.NetCDFVelocityReader}

restart_dict    = {'restart' : trajectory.RestartParser}

import topology
topology_dict  = {'amber'  : topology.TopologyParser}

converter_dict = {}
available_pottypes = ['amber-base', 'amber',  'amber94', 'amber96', 'amber99',
                      'amber99SB', 'amber03', 'amber10', 'amber12SB' ]
converter_dict = { ('amber', pot_type):topology.Format2AmberBaseConverter
                   for pot_type in available_pottypes }

coordinate_writer_dict = {
        'ascii' : trajectory.AsciiCoordinateWriter,
        'netcdf': trajectory.NetCDFCoordinateWriter}

velocity_writer_dict = {
        'ascii' : trajectory.AsciiVelocityWriter,
        'netcdf': trajectory.NetCDFVelocityWriter}
