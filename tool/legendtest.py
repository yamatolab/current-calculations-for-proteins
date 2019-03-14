import pygraphviz as pgv

filename = 'hoge.pdf'

dot = '''
digraph {
  rankdir=LR
  node [shape=plaintext]
  graph [splines=ortho]
  subgraph cluster_01 { 
    label = "Legend";
    key [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
      <tr><td align="right" port="i1">item 1</td></tr>
      <tr><td align="right" port="i2">item 2</td></tr>
      <tr><td align="right" port="i3">item 3</td></tr>
      <tr><td align="right" port="i4">item 4</td></tr>
      <tr><td align="right" port="i5">item 5</td></tr>
      </table>>]
    key2 [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0">
      <tr><td port="i1" bgcolor='greenyellow'>&nbsp;</td></tr>
      <tr><td port="i2">&nbsp;</td></tr>
      <tr><td port="i3">&nbsp;</td></tr>
      <tr><td port="i4">&nbsp;</td></tr>
      <tr><td port="i5">&nbsp;</td></tr>
      </table>>]
    key:i1:e -> key2:i1:w [color=red]
    key:i2:e -> key2:i2:w [color=gray]
    key:i3:e -> key2:i3:w [color=peachpuff3]
    key:i4:e -> key2:i4:w [color=turquoise4, style=dotted]
    key:i5:e -> key2:i5:w [color=red, style=dotted]
  }
}
'''

graph = pgv.AGraph(dot)
for n in graph.iteredges():
    print(n)
graph.add_node(1)
graph.add_node(2)
graph.add_edge(1, 2, penwidth=5)
graph.add_edge(1, 2, penwidth=1)
# print(graph.has_edge(2,1))
# print(graph.get_edge(1,2).attr['penwidth'])

graph.layout(prog='dot')
graph.draw(filename)

# a = {}
# a[(5,4)] = 5
