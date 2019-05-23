
import trajectory
coordinate_dict = {'presto-bin' : trajectory.CoordinateParser}
velocity_dict   = {'presto-bin' : trajectory.VelocityParser}
restart_dict    = {'presto-bin' : trajectory.RestartParser}

import topology
topology_dict  = {'presto'  : topology.TopologyParser}
converter_dict = {('presto', 'amber99') : topology.Format2Amber99Converter}

