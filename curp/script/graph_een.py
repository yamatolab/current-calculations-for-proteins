#! /usr/bin/env python
from __future__ import print_function
import os, sys
import itertools as it
import numpy
import pygraphviz as pgv

# FONTNAME = 'Helvetica Neue Bold'
FONTNAME = "Arial Bold"

class FormatWrongError(Exception): pass


def graph_een(data_fns, targets=['1-'], force_nodes=[], close_pairs=[],
              threshold=None, tpl_threshold=None, use_decrease=False,
              use_1letter=False, title='',
              ratio=None, alpha=False, graph_size=None, direction='TB',
              hide_isonode=False, line_values=[0.015, 0.008,0.003],
              line_colors=['red', 'blue', 'green'],
              line_thicks=[4.0, 4.0, 2.5], line_weights=[5.0, 3.0, 1.0],
              cluster_fn='', node_fn='', fig_fn=None, **kwds):

    if tpl_threshold is None:
        tpl_threshold = threshold

    # prepare
    parsers = [ EnergyConductivityParser(fp) for fp in data_fns ]

    # determine threshold value
    if threshold is None:
        threshold = line_values[-1]
    else:
        if threshold >= line_values[-1]:
            threshold = threshold
        else:
            threshold = line_values[-1]

    # communication chart
    cc  = CommunicationChart(
        title       = title if title else parsers[-1].get_title(),
        direction   = direction,
        threshold   = threshold,
        use_1letter = use_1letter,
        ratio       = ratio,
        alpha       = alpha,
        size        = graph_size,
        # font_scale  = 3.0 if opts.hide_node else 1.0,
        tpl_threshold = tpl_threshold,
        use_decrease = use_decrease)

    cc.set_line_thresholds(line_values,
            line_colors, line_thicks, line_weights)

    igrp_max = max([parser.get_max_id() for parser in parsers])
    igrp_min = min([parser.get_min_id() for parser in parsers])
    print('the range of id: {} - {}'.format(igrp_min, igrp_max))

    # donor and acceptor targets
    from curp.table.target import parse_target_atoms_line as parse_target_group
    if targets:
        targets = parse_target_group(targets, igrp_max)
        cc.set_targets(targets)

    if force_nodes:
        force_nodes = parse_target_group(force_nodes, igrp_max)
        cc.set_force_nodes(force_nodes)

    # if don_targets or opts.acc_targets:
        # don_targets = parse_target_group(opts.don_targets, igrp_max)
        # acc_targets = parse_target_group(opts.acc_targets, igrp_max)
        # cc.restrict_groups(don_targets, acc_targets)

    # traverse all parsers
    values = []
    for parser in parsers:
        vs = list(parser)
        values += vs
    vs_1 = vs # for parsers[-1]

    # create nodes and edges from all the data for defining graph topology
    # but not showing edge shapes.
    if hide_isonode:
        all_nodes = []
        for v in values:
            don, acc, cond = v[0], v[1], v[2]
            all_nodes.append(don)
            all_nodes.append(acc)
            cc.add_data(don, acc, cond)

        target_nodes = []
        for v in vs_1:
            don, acc, cond = v[0], v[1], v[2]
            if cond >= threshold:
                nid1, nname1, label1 = cc._parse_label(don)
                nid2, nname2, label2 = cc._parse_label(acc)
                target_nodes.append(nname1)
                target_nodes.append(nname2)

        other_nodes = set(all_nodes) - set(target_nodes)
        for n in other_nodes:
            cc.hide_node(n)

        # if opts.output_layout_fp:
            # cc.write_layout(opts.output_layout_fp)

        for v in vs_1:
            don, acc, cond = v[0], v[1], v[2]
            cc.show_edge(don, acc, cond)

        #BEGUG
        # print( " ".join(str(rid) for rid in sorted(set(target_nodes)) ))

    else:
        # create nodes and edges from all the data for defining graph topology
        # but not showing edge shapes.
        for v in values:
            don, acc, cond = v[0], v[1], v[2]
            cc.add_data(don, acc, cond)

        # if opts.output_layout_fp:
            # cc.write_layout(opts.output_layout_fp)

        # draw edges for parser[-1] (data_fps[-1])
        for v in vs_1:
            don, acc, cond = v[0], v[1], v[2]
            cc.show_edge(don, acc, cond)

    # apply node style
    if node_fn:

        node_parser = NodeStyleParser(node_fn, (igrp_min, igrp_max))
        for secname in node_parser.get_secnames():
            nodes = node_parser.get_vals(secname)
            attr_dict = node_parser.get_attr_dict(secname)
            cc.apply_node_attributes(rids=nodes, nnames=[], **attr_dict)

    # apply node pairs to close in together
    if close_pairs:

        for word in close_pairs:
            rid1, rid2 = word.split(':')
            cc.apply_close_nodes(rid1, rid2)

    # read clusters file and apply clusters to graph if clusters file exists.
    if cluster_fn:

        cls_parser = ClusterParser(cluster_fn, (igrp_min, igrp_max))
        for secname in cls_parser.get_secnames():
            rids = cls_parser.get_vals(secname)
            cc.create_cluster(label=secname, rids=rids)

            attr_dict = cls_parser.get_attr_dict(secname)
            cc.apply_cluster_attributes(secname, **attr_dict)

        # g.graph_attr.update(color='pink')

    cc.draw(fig_fn)



class EnergyConductivityParser:

    def __init__(self, filename):

        self.parse(filename)

    def __iter__(self):
        return iter(self._data_iter)

    def parse(self, filename):
        self._title, self._labels = self.parse_header_lines(filename)
        self._data_iter = list(self.gen_parse_data(filename))

    def _open(self, filename):
        if filename.endswith('.gz'):
            fd = gzip.open(filename, 'rb')
        else:
            fd = open(filename, 'rb')

        return fd

    def parse_header_lines(self, filename):
        """Parse header line."""

        fd = self._open(filename)

        title = ''
        labels = []

        lines = ( line for line in fd if line.startswith('#') )

        for line, iline in it.izip(lines, range(100)):

            if line.startswith('#title'):
                title = line.replace('#title', '')

            elif line.startswith('#label'):
                labels = line.split()[1:]

            else:
                pass

        fd.close()

        return title, labels

    def gen_parse_data(self, filename):

        # get lines
        fd = self._open(filename)
        lines = ( line for line in fd if not line.startswith('#') )

        for line in lines:
            cols = line.split()
            don, acc = cols[0], cols[1]
            cond = float(cols[2])
            others = cols[3:]

            yield don, acc, cond, others

        fd.close()

    def get_max_id(self):

        maxid = 0
        for values in self._data_iter:
            acc_res = values[1]
            rid = int(acc_res.split('_')[0])
            maxid = max(maxid, rid)

        return maxid

    def get_min_id(self):

        minid = 99999
        for values in self._data_iter:
            don_res = values[0]
            rid = int(don_res.split('_')[0])
            minid = min(minid, rid)

        return minid

    def get_title(self):
        return self._title



class CommunicationChart:

    _res_3_to_1 = dict(
            ALA='A', ARG='R', ASN='N', ASP='D', CYS='C', GLN='Q', GLU='E',
            GLY='G', HIS='H', ILE='I', LEU='L', LYS='K', MET='M', PHE='F',
            PRO='P', SER='S', THR='T', TRP='W', TYR='Y', VAL='V', HIP='H',
            HID='H', HIE='H', ASH='D', GLH='E', CYX='C',
            )
    __res_3_to_1_main = {k+'M':v+'M' for k,v in _res_3_to_1.items()}
    __res_3_to_1_side = {k+'S':v+'S' for k,v in _res_3_to_1.items()}

    _res_3_to_1.update(__res_3_to_1_main)
    _res_3_to_1.update(__res_3_to_1_side)

    def __init__(self, title,direction='TB', ratio=None,
            use_1letter=False, alpha=None, size="", **options):
        # self.__graph = pgv.AGraph(label=title, rankdir=direction)

        # preprocess the graph options
        graph_opts = {}

        if size:
            graph_opts["ratio"] = "auto"
            graph_opts["size"]  = size
        else:
            graph_opts["ratio"] = ratio if ratio else "auto"

        if alpha is not None:
            graph_opts["bgcolor"] = "#000000{:02}".format(alpha)

        self.__graph = pgv.AGraph(label=title, rankdir=direction, **graph_opts)

        # self._ec_1 = self._ec_1_with_alpha
        # self._ec_2 = self._ec_2_with_alpha
        # self._ec_3 = self._ec_3_with_alpha
        # self._ec_4 = self._ec_4_with_alpha

        self.__opts = options

        self._rid_to_nname = {}
        self.__label_to_subgraph = {}
        self._subcount = 0
        self.__use_1letter = use_1letter
        # ratio = ratio if ratio else 'auto'

        self.__targets = []
        self.__force_nodes = []
        self.__don_targets = []
        self.__acc_targets = []
        self._edge_to_cond = {}

        # set default attributes for graph
        self.__graph.graph_attr.update(
                fontsize = '14'     ,
                fontname = FONTNAME,
                splines  = 'polyline' ,
                ratio    = ratio,
                # splines  = 'ortho'    ,
                # size = "100x500",
        )

        # set default attributes for nodes
        self.__graph.node_attr.update(
                shape     = 'box'   ,
                style     = 'rounded,filled' ,
                fontsize  = '14'  ,
                fontname  = FONTNAME,
                fixedsize = True   ,
                # height    = 0.02,
                # width     = 0.02,
                fillcolor = 'grey'  ,
                fontcolor = 'black' ,
                color     = 'black'  ,
        )

        # set default attributes for edges
        self.__graph.edge_attr.update(
                fontname = FONTNAME,
        )

    def set_line_thresholds(self, values, colors, thicks, weights):

        self._ec_thresholds = []
        for v, c, t, w in zip(values, colors, thicks, weights):
            self._ec_thresholds.append( (float(v), str(c), str(t), str(w)) )

        # for transparency
        # self._ec_thresholds.append( (0.015,   '#FF000080'   , '8.0' , '5' ) )
        # self._ec_thresholds.append( (0.008,   '#0000FF80'  , '4.5' , '3' ) )
        # self._ec_thresholds.append( (0.003,  '#00800080' ,   '2.0' , '1' ) )
        # self._ec_thresholds.append( (0.001,  '#00000080' , '0.5' , '1' ))

    def apply_node_attributes(self, rids=[], nnames=[], **attr_dict):
        """Apply givin attrifutes to nodes."""

        nnames_all = [ nname for nname in self.__graph.nodes() ]
        rids_all = self._rid_to_nname.keys()

        for rid in rids:
            if rid not in rids_all: continue

            nname = self._rid_to_nname[rid]
            node = self.__graph.get_node(nname)
            node.attr.update(**attr_dict)

        for nname in nnames:
            if nname not in nnames_all: continue

            node = self.__graph.get_node(nname)
            node.attr.update(**attr_dict)

    def apply_ec_attributes(self, **attributes):
        self.__graph.edge_attr.update( **attributes )

    def _parse_label(self, line):
        rid, rname = line.split('_')

        if self.__use_1letter:
            rname = self._res_3_to_1.get(rname.upper(), rname)

        if self._is_special(line):
            nname = self._get_spetial_nname(rname)
            label = self._get_spetial_label(rname)
        elif rname[-1].isdigit():
            nname = '{}_{}'.format(rid,rname)
            label = '{}'.format(rname)
        else:
            nname = '{}_{}'.format(rid,rname)
            label = '{}{}'.format(rname.title(), int(rid))

        return int(rid), nname, label

    def _get_spetial_nname(self, rname):
        fst = rname.find('(')
        lst = rname.find(')')
        return rname[fst+1:lst]

    def _get_spetial_label(self, rname):
        fst = rname.find('(')
        lst = rname.find(')')
        return rname[:fst]+rname[fst+1:lst]+rname[lst+1:]
        # return rname.replace('(','').replace(')','')

    def _get_edge_attrs(self, *values):

        cond = values[0]
        tau = values[1] if len(values)>=2 else ''

        for ec_thresh_info in self._ec_thresholds:
            if cond >= ec_thresh_info[0]:
                color, width, weight = ec_thresh_info[1:]
                break

        # if cond >= self._ec_1[0]:
            # color, width, weight = self._ec_1[1:]
        # elif  cond >= self._ec_2[0]:
            # color, width, weight = self._ec_2[1:]
        # elif cond >= self._ec_3[0]:
            # color, width, weight = self._ec_3[1:]
        # else:
            # color, width, weight = self._ec_4[1:]

        return color, width, weight, tau

    def add_data(self, don, acc, cond, *others):

        rid1, nname1, label1 = self._parse_label( don )
        rid2, nname2, label2 = self._parse_label( acc )

        # if it isn't target groups then dont't add nodes
        is_forced1 = nname1 in self.__force_nodes
        is_forced2 = nname2 in self.__force_nodes

        if not (is_forced1 or is_forced2):
            if not (nname1 in self.__targets or nname2 in self.__targets):
                return
            else:
                if nname1 not in self.__don_targets: return
                if nname2 not in self.__acc_targets: return

        if cond >= self.__opts['tpl_threshold']:
            self.__graph.add_node( nname1, label=label1, pin=True)
            self.__graph.add_node( nname2, label=label2, pin=True)
            self._rid_to_nname[rid1] = nname1
            self._rid_to_nname[rid2] = nname2
            self.add_edge(nname1, nname2, cond)

        elif self.__opts['use_decrease'] \
                and cond <= -self.__opts['tpl_threshold']:
            self.__graph.add_node( nname1, label=label1, pin=True)
            self.__graph.add_node( nname2, label=label2, pin=True)
            self._rid_to_nname[rid1] = nname1
            self._rid_to_nname[rid2] = nname2
            self.add_edge(nname1, nname2, -cond)

        elif is_forced1:
            self.__graph.add_node( nname1, label=label1, pin=True )
            self._rid_to_nname[rid1] = nname1

        elif is_forced2:
            self.__graph.add_node( nname2, label=label2, pin=True )
            self._rid_to_nname[rid2] = nname2

        else:
            pass

    # def remove_node(self, node):
        # nid, nname, label = self._parse_label( node )
        # if self.__graph.has_node(nname):
            # print('ok', nname)
            # self.__graph.remove_node( nname )

    def _is_special(self, label):
        fst = label.find('(')
        lst = label.find(')')

        if fst>=0 and lst>=0:
            if fst<lst: return True
            else: raise FormatWrongError(label)
        elif fst<0 and lst<0:
            return False
        else:
            raise FormatWrongError(label)

    def hide_node(self, node):
        rid, nname, label = self._parse_label( node )
        if self.__graph.has_node(nname):
            node = self.__graph.get_node( nname )
            node.attr['style'] = 'invis'

    def remove_edges_from(self, nodes):
        node_names = []
        for n in nodes:
            rid, nname, label = self._parse_label( n )
            node_names.append(nname)

        self.__graph.remove_edges_from( [node_names] )

    def add_edge(self, nname1, nname2, cond, *others):

        cond_old = self._edge_to_cond.get((nname1, nname2))
        cond_old_rev = self._edge_to_cond.get((nname2, nname1))

        if cond_old:
            if cond > cond_old:
                self._edge_to_cond[key] = cond
                color, width, weight, edge_label = self._get_edge_attrs(cond)
                self.__graph.add_edge( nname1, nname2,
                        weight=weight, style='invis' )

        elif cond_old_rev:
            if cond > cond_old_rev:
                self._edge_to_cond[key] = cond
                color, width, weight, edge_label = self._get_edge_attrs(cond)
                self.__graph.add_edge( nname2, nname1,
                        weight=weight, style='invis' )

        else:
            self._edge_to_cond[key] = cond
            color, width, weight, edge_label = self._get_edge_attrs(cond)
            self.__graph.add_edge( key[0], key[1],
                    weight=weight, style='invis' )

        # if cond_old:
            # if cond > cond_old:
                # self._edge_to_cond[key] = cond
                # color, width, weight, edge_label = self._get_edge_attrs(cond)
                # self.__graph.add_edge( key[0], key[1],
                        # weight=weight, style='invis' )

        # else:
            # self._edge_to_cond[key] = cond
            # color, width, weight, edge_label = self._get_edge_attrs(cond)
            # self.__graph.add_edge( key[0], key[1],
                    # weight=weight, style='invis' )

    def apply_force(self):
        # if it isn't target groups then dont't add nodes
        is_forced1 = nname1 in self.__force_nodes
        is_forced2 = nname2 in self.__force_nodes

        if is_forced1:
            self.__graph.add_node( nname1, label=label1 )
            self._rid_to_nname[rid1] = nname1

        if is_forced2:
            self.__graph.add_node( nname2, label=label2 )
            self._rid_to_nname[rid2] = nname2

        # if not (is_forced1 or is_forced2):
            # if not (nname1 in self.__targets or nname2 in self.__targets):
                # return
            # else:
                # if nname1 not in self.__don_targets: return
                # if nname2 not in self.__acc_targets: return

    def show_edge(self, don, acc, cond, *others):

        thres = self.__opts['threshold']

        if -thres < cond < thres: return
        if cond < 0.0 and not self.__opts["use_decrease"]: return

        rid1, nname1, label1 = self._parse_label( don )
        rid2, nname2, label2 = self._parse_label( acc )

        if not self.__graph.has_edge(nname1, nname2):
            if not self.__graph.has_node(nname1):
                self.__graph.add_node(nname1, label=label1)
                self._rid_to_nname[rid1] = nname1
            if not self.__graph.has_node(nname2):
                self.__graph.add_node(nname2, label=label2)
                self._rid_to_nname[rid2] = nname2

            self.__graph.add_edge(nname1, nname2)

        edge = self.__graph.get_edge(nname1, nname2)

        node1 = self.__graph.get_node(nname1)
        node2 = self.__graph.get_node(nname2)
        node1.attr.update(self.__graph.node_attr)
        node2.attr.update(self.__graph.node_attr)
        node1.attr['label'] = label1
        node2.attr['label'] = label2

        # fscale = self.__opts['font_scale']
        # node1.attr['fontsize'] = fscale *

        if cond > 0.0:
            color, width, weight, edge_label = self._get_edge_attrs(cond)
            edge.attr['style'] = ''
            edge.attr['color'] = color
            edge.attr['width'] = width
            edge.attr['label'] = edge_label
            edge.attr['penwidth'] = width
        else:
            color, width, weight, edge_label = self._get_edge_attrs(-cond)
            edge.attr['style'] = 'dashed'
            edge.attr['color'] = color
            edge.attr['width'] = width
            edge.attr['label'] = edge_label
            edge.attr['penwidth'] = width

    def restrict_groups(self, don_targets, acc_targets):
        self.__don_targets = don_targets
        self.__acc_targets = acc_targets

    def set_targets(self, targets):
        self.__targets = targets

    def set_force_nodes(self, force_nodes):
        self.__force_nodes = force_nodes

    def create_cluster(self, label, rids, **attributes):
        """Create subgraph cluster by label and nodes."""

        # default attributes
        default_attributes = dict(
                fontsize = 24,
                style = 'rounded,filled',
                color = 'lightgrey', )

        default_attributes.update( **attributes )

        # if node exists in graph
        # united_nnames = list( set(self._rid_to_nname.keys()) & set(rids) )
        united_nnames = [ self._rid_to_nname[rid]
                for rid in rids if self._rid_to_nname.get(rid) ]

        # create new graph
        self._subcount += 1
        name = 'cluster' + str(self._subcount)
        # name = 'a' + str(self._subcount)
        g = self.__graph.subgraph(nbunch=united_nnames, # rank='same',
                name=name, label=label, **default_attributes)
        self.__label_to_subgraph[label] = g
        # g.graph_attr['rank'] = 'same'
        return g

    def apply_cluster_attributes(self, secname, **attributes):
        g = self.__label_to_subgraph[secname]
        g.graph_attr.update( **attributes )

    def write_layout(self, fp):
        self.__graph.layout(prog='dot')
        self.__graph.write(fp)

    def draw(self, filename):
        # print('create legend')
        # self.create_legend()
        print('do layout')
        # self.__graph.layout(prog='dot')
        # self.__graph.write('tmp.dot')
        # self.__graph.write('test1.dot')
        # self.__graph.write('test2.dot')
        print('do draw')
        self.__graph.draw(filename, prog='dot')

    def apply_close_nodes(self, rid1, rid2, weight=3):

        nname1 = self._rid_to_nname.get(rid1)
        nname2 = self._rid_to_nname.get(rid2)

        if nname1 is None: return
        if nname2 is None: return

        if self.__graph.has_edge(nname1, nname2):
            edge = self.__graph.get_edge(nname1, nname2)
            edge.attr['weight'] = weight
        else:
            self.__graph.add_edge( nname1, nname2,
                    weight=weight, style='invis' )

    def get_nodes(self):
        return sorted(self._rid_to_nname.values())

    def create_legend(self):
        """Create a legend using the function of the subgraph."""

        # default attributes
        default_attributes = dict(
                fontsize = 48,)
                # style = 'rounded,filled',)
                # color = 'grey', fillcolor='white')

        # default_attributes.update( **attributes )
        # label = """
    # <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
     # <TR>
      # <TD COLSPAN="2"><B>Legend</B></TD>
     # </TR>
     # <TR>
      # <TD>Foo</TD>
      # <TD><FONT COLOR="red">Foo</FONT></TD>
     # </TR>
     # <TR>
      # <TD>Bar</TD>
      # <TD BGCOLOR="RED"></TD>
     # </TR>
     # <TR>
      # <TD>Baz</TD>
      # <TD BGCOLOR="BLUE"></TD>
     # </TR>
     # <TR>
      # <TD>Test</TD>
      # <TD><IMG src="so.png" SCALE="False" /></TD>
     # </TR>
     # <TR>
      # <TD>Test</TD>
      # <TD CELLPADDING="4">
       # <TABLE BORDER="1" CELLBORDER="0" CELLSPACING="0" CELLPADDING="0">
        # <TR>
         # <TD BGCOLOR="Yellow"></TD>
        # </TR>
       # </TABLE>
      # </TD>
     # </TR>
    # </TABLE>

        # """

        # label = '''<
        # <font color="red">aaa</font> <font color="blue">bbb</font>
        # >'''
        label_key1 = '''<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
        <tr><td align="right" port="i1">item 1</td></tr>
        <tr><td align="right" port="i2">item 2</td></tr>
        <tr><td align="right" port="i3">item 3</td></tr>
        <tr><td align="right" port="i4">item 4</td></tr>
        </table>>'''

        label_key2 = '''<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
      <tr><td port="i1">&nbsp;</td></tr>
      <tr><td port="i2">&nbsp;</td></tr>
      <tr><td port="i3">&nbsp;</td></tr>
      <tr><td port="i4">&nbsp;</td></tr>
      </table>>'''

        self.__graph.add_node('key1', label=label_key1, color='white',
                fillcolor='white', fontsize=32)
        self.__graph.add_node('key2', label=label_key2)
        # self.__graph.add_node('key1:i1:e')
        # self.__graph.add_node('key2:i1:w')

        for n in self.__graph.iternodes():
            print(n)
        # self.__graph.add_node('key1:1', label='hoge1')
        # self.__graph.add_node('key2:1', label='', fillcolor='white')

        # self.__graph.add_node('key1:2', label='hoge2')
        # self.__graph.add_node('key2:2', label='', fillcolor='white')
        # print(self.__graph.get_node('key:i1'))


        # self.__graph.add_edge('key1', 'key2', style='dashed',
                # **default_attributes)
        # self.__graph.add_edge('key1', 'key2', style='dashed', color='blue',
                # penwidth=10, **default_attributes)

        self.__graph.add_edge('key1:2', 'key2:2', penwidth=10,
                **default_attributes)


        # if node exists in graph
        # label = 'hoge'
        # self.__graph.add_node( 'legend', label=label, **default_attributes )

        # create new graph
        # name = 'legend' + str(self._subcount)
        name = 'legend'
        # # name = 'a' + str(self._subcount)
        # print(label)
        # nodes = ['key1:1', 'key1:2', 'key2:1', 'key2:2']
        nodes = ['key1', 'key2']
        g = self.__graph.subgraph(nbunch=nodes, rank='sink',
                name=name, label='Legend', rankdir='LR')# **default_attributes)
        # self.__label_to_subgraph[name] = g
        # # g.graph_attr['rank'] = 'sink'
        # return g


from curp.table.ini_parser import IniParser
class ClusterParser(IniParser):

    def __init__(self, filename, igrp_range, *args, **kwds):
        self._igrp_range = igrp_range
        IniParser.__init__(self, filename, *args, **kwds)

    def _gen_col(self, lines):
        groups = []
        for line in lines:
            if '=' in line: continue
            for col in line.split():

                if '-' in col:
                    groups.extend( self._parse_mask_with_range(col, groups) )

                else:
                    groups.append( int(col) )

        return groups

    def _parse_mask_with_range(self, mask, groups):

        cols = mask.split('-')

        if cols[0] != '' and cols[1] != '':
            igrp_beg, igrp_end = int(cols[0]), int(cols[1])

        elif cols[0] == '' and cols[1] != '': # -45
            igrp_end = int(cols[1])

            if len(groups) == 0:
                igrp_beg = self._igrp_range[0]
            else:
                igrp_beg = groups[-1] + 1

        elif cols[0] != '' and cols[1] == '': # 45-
            igrp_beg = int(cols[0])
            igrp_end = self._igrp_range[1]

        else:
            pass

        return range(igrp_beg, igrp_end+1)

# define NodeStyleParser from ClusterParser
NodeStyleParser = ClusterParser


if __name__ == '__main__':
    from console import arg_graph_een, exec_command

    parser = arg_graph_een()
    exec_command(parser)
