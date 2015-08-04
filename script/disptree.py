from ete2 import PhyloNode, TreeStyle
import sys
if(len(sys.argv)<2):
    sys.exit("Usage : python disptree.py input output")
else:
    inputfile, output = sys.argv[1], sys.argv[2]
    t = PhyloNode(inputfile)
    ts = TreeStyle()
    ts.show_branch_length = True
    ts.show_branch_support = True
    if not output.endswith('.pdf'):
        output+='.pdf'
    t.render(output, tree_style=ts)
