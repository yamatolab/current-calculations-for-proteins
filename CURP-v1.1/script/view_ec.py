from pymol import cmd
import sys, os

ONE_LETTER = {
        'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',
        'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',
        'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',
        'GLY':'G', 'PRO':'P', 'CYS':'C',
        'HIE':'H', 'HID':'H', 'HIP':'H', 'ASH':'D', 'GLH':'E',
        'WAT':'w', 'HOH':'w',
        }

def main():

    args = get_args()

    ecs = list(parse_ec(args.ec_fps[0]))

    inc_ecs, dec_ecs = [], []
    for data in ecs:
        ec = float(data['values'][0])
        if ec >= args.threshold:
            inc_ecs += [(data, ec)]
        elif ec <= - args.threshold:
            dec_ecs += [(data, ec)]
        else:
            pass
    
    if args.use_increase:
        for d, ec in inc_ecs:
            display_pair(d, ec, args)

    if args.use_decrease:
        for d, ec in dec_ecs:
            display_pair(d, ec, args)

    if args.output_fp:
        cmd.ray(args.resolution)
        cmd.png(args.output_fp)
        cmd.quit()

    if args.movie_style:

        # set resolution
        x, y = cmd.get_viewport(output=1)
        cmd.viewport(args.resolution, args.resolution*y/x)

        # perform ray trace in write the movie
        cmd.set("ray_trace_frames", 1)

        if args.movie_style == 'y-rotation':
            cmd.mset(1, 120)
            cmd.movie.roll(1,120,1,axis='y')
            cmd.mpng("_movie/")
            cmd.quit()

        elif args.movie_style == 'x-rotation':
            cmd.mset(1, 120)
            cmd.movie.roll(1,120,1,axis='x')
            cmd.mpng("_movie/")
            cmd.quit()

        elif args.movie_style == 'y-rock':
            cmd.mset(1, 120)
            cmd.movie.rock(1, 120, 40, axis='y')
            cmd.mpng("_movie/")
            cmd.quit()

        elif args.movie_style == 'x-rock':
            cmd.mset(1, 120)
            cmd.movie.rock(1, 120, 40, axis='x')
            cmd.mpng("_movie/")
            cmd.quit()

        elif args.movie_style == 'nutate':
            cmd.mset(1, 120)
            cmd.movie.nutate(1, 120, 40)
            cmd.mpng("_movie/")
            cmd.quit()

        else:
            pass


def get_args():
    """Parse and get the command line arguments."""
    import argparse
    parser = argparse.ArgumentParser(
            description="Show 3D structure with EEN.")

    parser.add_argument(
            dest='ec_fps', nargs='+', metavar='EC_FP',
            help="specify the ec file paths."
            )

    # parser.add_argument(
            # '-r', "--important-residues", dest='residues',
            # metavar=('res1','res2'), default=[], type=int, nargs='+',
            # help="specify residues that you want to highlight.",
            # )

    parser.add_argument(
            "-L", "--without-residue-labels", dest="enable_labels",
            action="store_false",
            help="enable the labels of residues."
            )

    parser.add_argument(
            "-I", "--not-dispaly-increase", dest="use_increase",
            action="store_false",
            help="don't display increasing pairs."
            )

    parser.add_argument(
            "-d", "--dispaly-decrease", dest="use_decrease",
            action="store_true",
            help="display dncreasing pairs."
            )

    parser.add_argument(
            "-t", "--threshold", dest="threshold",
            default=0.003, type=float,
            help="specify the threshold value that show on 3D structure.",
            )

    parser.add_argument(
            "-n", "--node-color", dest="node_color",
            default='grey', type=str,
            help="specify the color of the residues."
            )

    parser.add_argument(
            "-o", "--output-filepath", dest="output_fp",
            default="",
            help="specify a file path to write figure."
            )

    parser.add_argument(
            "-r", "--resolution", dest="resolution",
            default=820, type=int,
            help="specify the resolution of a file path to write figure."
            )

    # parser.add_argument(
            # '-m', "--movie-filename", dest="movie_fp",
            # default="",
            # help="specify a file path to write a movie."
            # )

    parser.add_argument(
            '-ms', "--movie-style", dest="movie_style",
            default="y-rotaiton",
            choices=["y-rotation", "x-rotation", "x-rock", "y-rock", "nutate"],
            help="specify the style of the movie."
            )

    args = parser.parse_args()
    return args

def parse_ec(ec_fp):

    file = open(ec_fp, 'r')

    lines = (line.strip() for line in file
            if line
            if not line.startswith('#') )

    for line in lines:

        cols = line.split()
        res1, res2, values = cols[0], cols[1], cols[2:]

        rid1, rname1 = res1.split('_')
        rid2, rname2 = res2.split('_')

        yield dict(
                rid1=int(rid1), rname1=rname1,
                rid2=int(rid2), rname2=rname2,
                values = values )

    file.close()

def highlight_important_residues(filename, residue_shift=0):

    colors = ['red', 'orange', 'yellow', 'green']

    file = open(filename, 'b')

    # remove the only space line and the comment line
    line_iter = ( line.strip() for line in file
            if not line.isspace() if not line.startswith('#') )

    for line in line_iter:
        cols = line.split()
        rid, irank = int(cols[0]), int(cols[1])
        sel = str(rid-residue_shift) + '/CA'

        # show spheres
        cmd.show('spheres', sel)
        cmd.label(sel)

        # coloring to spheres
        cmd.set('sphere_color', colors[irank-1], sel)
        cmd.set('sphere_scale', 0.6, sel)

    file.close()


def get_selection(rid, rname):

    if rname == 'WAT':
        sel = '(resi ' + str(rid) + ' and name O)'
    elif rname == 'CAU':
        sel = '(resn CAU and name O2)'
    else:
        sel = str(rid) + '/CA'

    return sel

def display_pair(data, ec, args):

    is_inc = ec > 0
    ec = abs(ec)

    sel1 = get_selection(data['rid1'], data['rname1'])
    sel2 = get_selection(data['rid2'], data['rname2'])

    # connect pair by ec value
    rn1 = ONE_LETTER.get(data['rname1'])
    rn2 = ONE_LETTER.get(data['rname2'])
    rid1 = str(data['rid1'])
    rid2 = str(data['rid2'])
    name1 = rn1+rid1 if rn1 else rid1
    name2 = rn2+rid2 if rn2 else rid2
    name  = "{}-{}".format(name1, name2)

    cmd.distance(name, selection1=sel1, selection2=sel2)
    cmd.hide("label", name)

    if 0.003 <= ec < 0.008:
        cmd.set('dash_width', 4, name)
        # cmd.set('dash_color', 'orange', name)
        cmd.set('dash_color', 'green', name)
    elif 0.008 <= ec < 0.015:
        cmd.set('dash_width', 4, name)
        cmd.set('dash_color', 'blue', name)
        # cmd.set('dash_color', 'white', name)
    elif 0.015 <= ec:
        cmd.set('dash_width', 6, name)
        cmd.set('dash_color', 'red', name)
    else:
        pass

    cmd.set('dash_gap', 0.8, name)
    if args.use_increase and args.use_decrease:
        if is_inc: cmd.set('dash_gap', 0.0, name)

    # show residues consisting of the pair 
    # cmd.show('sticks', sel1.replace('/CA','/'))
    # cmd.show('sticks', sel2.replace('/CA','/'))
    # cmd.color('blue', sel1.replace('/CA','/'))

    cmd.show('spheres', sel1)
    cmd.show('spheres', sel2)
    cmd.set('sphere_scale', 0.4, sel1)
    cmd.set('sphere_scale', 0.4, sel2)
    cmd.color(args.node_color, sel1)
    cmd.color(args.node_color, sel2)

    # label the pair
    if args.enable_labels:
        cmd.label(sel1, '"{}"'.format(name1))
        cmd.label(sel2, '"{}"'.format(name2))
        cmd.set("label_size", 20)
        # cmd.set("label_position", (-1.0,-1.0,-1.0))

    else:
        cmd.hide('labels')

if __name__ == '__main__':
    main()
