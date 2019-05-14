import calvolume as cv

def get_volume_calculator(topology, setting, target_atoms, gname_iatoms_pairs):
    """Return VolumeCalculator object by setting."""

    # get volume setting
    traj_parser = None
    vsetting = cv.VolumeSetting(topology, setting, traj_parser,
                                target_atoms, gname_iatoms_pairs)

    method = setting.volume.method
    if method == 'none':
        obj = cv.Volume1(vsetting)
    elif method == 'vdw':
        obj = cv.VDWVolumeCalculator(vsetting)
    # elif method == 'smve': # provided, by Srolovitz, Maeda, Vitek, Egami
    #     obj = cv.SMVEVolumeCalculator(vsetting)
    elif method == 'voronoi': # use voronoi polyhedron
        if vsetting.voronoi_solvation == 'none':
            obj = cv.VoronoiVolumeCalculator(vsetting)

        else: # use solvation
            obj = cv.VoronoiVolumeCalculatorWithSolvation(vsetting)

    elif method == 'outer': # use outer volume trajectory already calculated
        obj = cv.OuterVolumeFetcher(vsetting)
    else:
        obj = cv.Volume1(vsetting)

    return obj
