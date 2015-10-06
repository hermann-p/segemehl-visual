def coffmanGraham(g):
    print 'Assigning layers with Coffman-Graham algorithm'
    lg = {}
    d, verts = transitiveReduction(g)
    # find starting element
    done = 0
    for v in verts:
        if len(v.left) == 0:
            v.layer = 0
            v.processed = True
            done += 1
    while done < Node.N:
        right = 0
        for i in xrange(Node.N):
            good = True
            right = 0
            for j in xrange(i):
                if d[i][j] == 1:
                    if not verts[j].processed:
                        good = False
                        break
                    right = max(right, verts[j].layer)
            if good: # all incoming/left nodes positioned
                verts[i].layer = right + 1
                verts[i].processed = True
                done += 1
    lg = {}
    for el in verts:
        L = el.layer
        if not L in lg.keys():
            lg[L] = [el]
        else:
            lg[L].append(el)
    return lg, verts
    