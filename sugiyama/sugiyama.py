import Tkinter as tk
from Queue import Queue
from collections import defaultdict

class Node:
    N = 0
    def __init__(self, name='', left=[], right=[], dummy=False):
        self.name = name
        self.left = []
        self.right = []
        for n in left:
            self.connectLeft(n)
        for n in right:
            self.connectRight(n)
        self.layer = None
        self.y = 0
        self.processed = False
        self.dummy = dummy
	self.id = Node.N
	Node.N += 1
	
    def connectLeft(self, lt):
        self.left.append(lt)
        lt.right.append(self)

    def connectRight(self, rt):
        self.right.append(rt)
        rt.left.append(self)

    def markUnprocessed(self):
        self.processed = False
        for n in self.left:
            if n.processed: n.markUnprocessed()
        for n in self.right:
            if n.processed: n.markUnprocessed()

def initNodes(testcase=0):
    print 'Setting up graph...'
    global graph
    if testcase == 0:
	graph = Node('A')
	C = Node('C', left=[graph])
	B = Node('B', left=[graph, C])
	D = Node('D', left=[B,C])
	E = Node('E', left=[B])
	F = Node('F', left=[D, E])
	H = Node('H', left=[F])
	G = Node('G', left=[D, H])
    elif testcase == 1:
	graph = Node('A')
	B = Node('B', left=[graph])
	C = Node('C', left=[B])
	D1 = Node('D1', left=[C])
	D2 = Node('D2', left=[C])
    elif testcase == 2:
	graph = Node('D1')
	C = Node('C', right=[graph])
	D2 = Node('D2', left=[C])
	B = Node('C', right=[C])
	A = Node('A', right=[B])
    else:
	graph = Node('A')
	B = Node('B', right=[graph])
	C = Node('C', left=[graph])
	D = Node('D', right=[graph])
	E = Node('E', left=[C])
	F = Node('F', right=[graph, B, C])
    
def assignLayers(node, layeredGraph, L=0, process=True):
    print 'Assigning layers...'
    if process:
	node.layer = L
	node.processed = True
	if L not in layeredGraph.keys():
	    layeredGraph[L] = []
	layeredGraph[L].append(node)
	print '-- assigned node', node.name, 'to layer', L
    for n in node.left:
        if not n.processed: assignLayers(n, layeredGraph, L-1)
    for n in node.right:
        if not n.processed: assignLayers(n, layeredGraph, L+1)

def assignLayers2(g):
    lg = {}
    lg = defaultdict(lambda: [], lg)
    tr, verts = transitiveReduction(g)
    g.layer = 0
    g.processed = True
    assigned = 1
    lg[0].append(g)
    while assigned < Node.N:
	for i in xrange(Node.N):
	    if verts[i].processed:
		print '-- node #', i, ':', verts[i].name, 'assigned, searching links'
		for j in xrange(i, Node.N):
		    if not verts[j].processed and tr[i][j] == 1:
			verts[j].layer = verts[i].layer + 1
			verts[j].processed = True
			assigned += 1
			lg[verts[j].layer].append(verts[j])
			print '++ assigned node', verts[j].name, 'to layer', verts[j].layer, ';', Node.N - assigned, 'left'
    g.markUnprocessed()
    return lg

def transitiveReduction(g):
    print 'Calculating transitive reduction...'
    d = [ [0 for i in xrange(Node.N)] for i in xrange(Node.N)]
    v = [None for i in xrange(Node.N)]
    q = Queue()
    q.put(g)
    while not q.empty():
        node = q.get()
        node.processed = True
        v[node.id] = node
        for n in node.left:
            if not n.processed: q.put(n)
            d[node.id][n.id] = d[n.id][node.id] = 1
        for n in node.right:
            if not n.processed: q.put(n)
            d[n.id][node.id] = d[node.id][n.id] = 1
    g.markUnprocessed()
    for x in v:
        for y in v:
            for z in v:
                if d[x.id][y.id] == 1 and d[y.id][z.id] == 1:
                    d[z.id][x.id] = d[x.id][z.id] = 0
    return d, v

def barycenters(lg):
    print 'Calculating y positions with barycenter heuristics...'
    old = Node.N
    verts = [None for n in xrange(Node.N)]
    print '-- Real nodes:', old
    insertDummies(lg)
    print '-- Dummy nodes:', Node.N - old
    tr, _ = transitiveReduction(lg[0][0])
    for lin in tr:
	print lin
    for L in sorted(lg.keys()):
	print '-- processing layer', L, 'with', len(lg[L]), 'elements'
	for n in lg[L]:
	    if not n.dummy: verts[n.id] = n
	    n.y = 0
	    pred = 0.0
	    for lt in n.left:
		if tr[n.id][lt.id] == 1:
		    print '---- adding', lt.name, lt.id, 'value:', lt.y
		    n.y += lt.y
		    pred += 1.0
	    if pred > 0:
		print '---- dividing by', pred
		n.y = int( n.y / pred )
        print 'OFFSET:', len(lg[L]), len(lg[min(lg.keys())]), len(lg[min(lg.keys())]) - len(lg[L])
	offset = len(lg[min(lg.keys())]) - len(lg[L])
	for i in xrange( 1, len(lg[L]) ):
	    if lg[L][i].y == lg[L][i-1].y:
		offset += 2
	    lg[L][i].y += offset
	for i in xrange(len(lg[L])):
	    print '++', L, '#', i, 'y:', lg[L][i].y, 'name:', lg[L][i].name if lg[L][i].name else 'dummy'
    return verts

def paintNodes(window, g):
    print 'Displaying nodes...'
    sx = 40
    sy = 40
    for node in g:
	if not node or node.dummy: continue
        print node.name, 'x:', node.layer, 'y:', node.y
        layer = node.layer
        for rt in node.right:
	    if rt.dummy: continue
	    print '-- drawing:', node.name, str(node.layer)+'/'+str(node.y), '->', rt.name, str(rt.layer)+'/'+str(rt.y)
            window.create_line((layer + 2) * sx, node.y * sy,
                               (rt.layer + 2) * sx, rt.y * sy)
	window.create_oval((layer + 1.8) * sx, (node.y + 0.2) * sy,
			    (layer + 2.2) * sx, (node.y - 0.2) * sy, fill='red')
	window.create_text( ((layer+2) * sx, node.y * sy), justify=tk.CENTER, text=node.name)
            
def insertDummies(lg):
    print 'Filling spans with dummy nodes'
    for L in sorted(lg.keys()):
	for n in lg[L]: # all nodes in layer
	    for rt in n.right:
		if rt.layer - n.layer > 1:
		    print 'layerspan:', n.name, n.layer, '->', rt.layer, rt.name
		    i = n.layer + 1
		    prev = n
		    while i < rt.layer:
			print '-- adding dummy level', i
			new = Node(left=[prev], dummy = True)
			lg[i].append(new)
			new.layer = i
			prev = new
			i += 1
		    prev.connectRight(rt)
		    print 'connecting dummy', i-1, 'to', rt.layer

def clearDummies(lg):
    print 'Removing dummy nodes'
    for L in lg.keys():
	for n in lg[L]:
	    if n.dummy: n.unDummy()
	    
def centerGraph(lg):
    lt = lg[0][0].layer
    tp = lg[0][0].y
    for L in lg:
	for v in lg[L]:
	    if not v: continue
	    lt = min(v.layer, lt)
	    tp = min(v.y, tp)
    xoffs = 2 - lt
    yoffs = 2 - tp
    for L in lg:
	for vv in lg[L]:	
	    if not vv: continue
	    vv.layer += xoffs
	    vv.y += yoffs
    


graph = None
layeredGraph = {}
    
main = tk.Tk()
window = tk.Canvas(main, width=600, height=450, background='white')
window.pack()

initNodes(testcase=0)
#assignLayers(graph, layeredGraph)
layeredGraph = assignLayers2(graph)
#layeredGraph, verts = coffmanGraham(graph)
for i in xrange(1):
    verts = barycenters(layeredGraph)
#clearDummies(layeredGraph)
centerGraph(layeredGraph)
paintNodes(window, verts)

tk.mainloop()
