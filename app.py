#!/usr/bin/env python3
"""
Graph Theory Lab — Modern UI Edition
Logistics & Urban Planning

All algorithms implemented in pure Python.
Run:  python app.py
Open: http://localhost:8765
"""

import http.server
import json
import heapq
import os
import sys
import threading
import webbrowser
from collections import defaultdict, deque

PORT = 8765
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# ════════════════════════════════════════════════════════════
#  GRAPH UTILITIES
# ════════════════════════════════════════════════════════════

def _build_adj(vertices, edges, directed):
    v_ids = {v['id'] for v in vertices}
    adj = defaultdict(list)
    for e in edges:
        u, v = e['from'], e['to']
        if u in v_ids and v in v_ids:
            adj[u].append((v, e.get('weight', 1), e.get('id', '')))
            if not directed:
                adj[v].append((u, e.get('weight', 1), e.get('id', '')))
    return adj, v_ids


def _fmt(val):
    """Format a number for display — int if whole, else float."""
    if val is None or val == float('inf'):
        return '∞'
    if isinstance(val, float) and val == int(val):
        return str(int(val))
    return str(val)


# ════════════════════════════════════════════════════════════
#  1. BASIC ANALYSIS
# ════════════════════════════════════════════════════════════

def analyze_graph(data):
    vertices = data['vertices']
    edges    = data['edges']
    directed = data.get('directed', False)

    n = len(vertices)
    m = len(edges)
    v_ids   = [v['id'] for v in vertices]
    v_index = {vid: i for i, vid in enumerate(v_ids)}

    degrees     = {vid: 0 for vid in v_ids}
    in_degrees  = {vid: 0 for vid in v_ids}
    out_degrees = {vid: 0 for vid in v_ids}

    for e in edges:
        u, v = e['from'], e['to']
        if directed:
            out_degrees[u] += 1
            in_degrees[v]  += 1
        degrees[u] += 1
        degrees[v] += 1

    # Connectivity (undirected BFS)
    adj_un = defaultdict(set)
    for e in edges:
        adj_un[e['from']].add(e['to'])
        adj_un[e['to']].add(e['from'])

    visited    = set()
    components = []
    for vid in v_ids:
        if vid not in visited:
            comp  = []
            queue = deque([vid])
            visited.add(vid)
            while queue:
                cur = queue.popleft()
                comp.append(cur)
                for nb in adj_un[cur]:
                    if nb not in visited:
                        visited.add(nb)
                        queue.append(nb)
            components.append(comp)

    has_neg = any(e.get('weight', 1) < 0 for e in edges)

    # Matrices
    adj_mat = [[0]*n for _ in range(n)]
    w_mat   = [[0]*n for _ in range(n)]
    for e in edges:
        i = v_index.get(e['from'])
        j = v_index.get(e['to'])
        if i is not None and j is not None:
            adj_mat[i][j] = 1
            w_mat[i][j]   = e.get('weight', 1)
            if not directed:
                adj_mat[j][i] = 1
                w_mat[j][i]   = e.get('weight', 1)

    return {
        'vertexCount':       n,
        'edgeCount':         m,
        'degrees':           degrees,
        'inDegrees':         in_degrees if directed else None,
        'outDegrees':        out_degrees if directed else None,
        'sumDegrees':        sum(degrees.values()),
        'connected':         len(components) <= 1,
        'componentCount':    len(components),
        'components':        components,
        'hasNegativeWeight': has_neg,
        'adjacencyMatrix':   adj_mat,
        'weightMatrix':      w_mat,
        'vertexLabels':      v_ids,
        'directed':          directed,
    }


# ════════════════════════════════════════════════════════════
#  2. DIJKSTRA  —  Shortest Path
# ════════════════════════════════════════════════════════════

def dijkstra(data):
    vertices = data['vertices']
    edges    = data['edges']
    directed = data.get('directed', False)
    source   = data['source']
    target   = data.get('target')

    if any(e.get('weight', 1) < 0 for e in edges):
        return {'error': 'Dijkstra ไม่รองรับน้ำหนักติดลบ!', 'steps': []}

    adj, v_ids = _build_adj(vertices, edges, directed)

    INF  = float('inf')
    dist = {v: INF for v in v_ids}
    prev = {v: None for v in v_ids}
    pedge = {v: None for v in v_ids}
    dist[source] = 0
    visited = set()
    steps   = []

    steps.append({
        'type': 'init',
        'desc': f'เริ่มต้น: dist[{source}] = 0, จุดอื่นทั้งหมด = ∞',
        'distances': {k: (v if v != INF else None) for k, v in dist.items()},
        'visited': [], 'current': source,
        'hlVerts': [source], 'hlEdges': [], 'pathEdges': [],
    })

    pq = [(0, source)]
    while pq:
        d, u = heapq.heappop(pq)
        if u in visited:
            continue
        visited.add(u)

        steps.append({
            'type': 'visit',
            'desc': f'เยี่ยมจุด {u}  (dist = {_fmt(d)})',
            'distances': {k: (v if v != INF else None) for k, v in dist.items()},
            'visited': list(visited), 'current': u,
            'hlVerts': [u], 'hlEdges': [], 'pathEdges': [],
        })

        if target and u == target:
            break

        for v, w, eid in adj[u]:
            if v not in visited:
                nd = d + w
                if nd < dist[v]:
                    old = dist[v]
                    dist[v]  = nd
                    prev[v]  = u
                    pedge[v] = eid
                    heapq.heappush(pq, (nd, v))

                    steps.append({
                        'type': 'relax',
                        'desc': (f'Relax ({u} → {v}): dist[{v}] = '
                                 f'{_fmt(old)} → {_fmt(nd)}  (ผ่าน {u}, w={_fmt(w)})'),
                        'distances': {k: (vv if vv != INF else None)
                                      for k, vv in dist.items()},
                        'visited': list(visited), 'current': u,
                        'considering': v,
                        'hlVerts': [u, v], 'hlEdges': [eid], 'pathEdges': [],
                    })

    # Build path
    path = []
    path_eids = []
    if target and dist.get(target, INF) != INF:
        c = target
        while c is not None:
            path.append(c)
            if pedge[c]:
                path_eids.append(pedge[c])
            c = prev[c]
        path.reverse(); path_eids.reverse()

        steps.append({
            'type': 'result',
            'desc': f'เส้นทางสั้นที่สุด: {" → ".join(path)},  ระยะรวม = {_fmt(dist[target])}',
            'distances': {k: (v if v != INF else None) for k, v in dist.items()},
            'visited': list(visited), 'path': path,
            'hlVerts': path, 'hlEdges': [], 'pathEdges': path_eids,
            'totalDistance': dist[target],
        })
    elif target:
        steps.append({
            'type': 'result',
            'desc': f'ไม่มีเส้นทางจาก {source} ไปยัง {target}',
            'distances': {k: (v if v != INF else None) for k, v in dist.items()},
            'visited': list(visited), 'path': [],
            'hlVerts': [], 'hlEdges': [], 'pathEdges': [],
            'totalDistance': None,
        })
    else:
        steps.append({
            'type': 'result',
            'desc': f'Dijkstra จาก {source} เสร็จสมบูรณ์ — คำนวณระยะทางทั้งหมดแล้ว',
            'distances': {k: (v if v != INF else None) for k, v in dist.items()},
            'visited': list(visited), 'path': [],
            'hlVerts': list(visited), 'hlEdges': [], 'pathEdges': [],
            'totalDistance': None,
        })

    return {'steps': steps}


# ════════════════════════════════════════════════════════════
#  3. PRIM  —  Minimum Spanning Tree / Forest
# ════════════════════════════════════════════════════════════

def prim(data):
    vertices = data['vertices']
    edges    = data['edges']
    start    = data.get('start', vertices[0]['id'] if vertices else None)

    if not vertices:
        return {'steps': [], 'totalWeight': 0}

    # Always undirected for MST
    adj = defaultdict(list)
    v_ids = {v['id'] for v in vertices}
    for e in edges:
        u, v = e['from'], e['to']
        if u in v_ids and v in v_ids:
            adj[u].append((v, e.get('weight', 1), e.get('id', '')))
            adj[v].append((u, e.get('weight', 1), e.get('id', '')))

    steps = []
    mst_eids   = []
    total_w    = 0
    visited    = set()
    remaining  = set(v_ids)
    roots      = []

    def _prim_comp(sv):
        nonlocal total_w
        visited.add(sv); remaining.discard(sv); roots.append(sv)

        steps.append({
            'type': 'init',
            'desc': f'เริ่ม Prim จากจุด {sv}',
            'visited': list(visited), 'mstEdges': list(mst_eids),
            'hlVerts': [sv], 'hlEdges': [], 'totalWeight': total_w,
        })

        heap = []
        for to_v, w, eid in adj[sv]:
            heapq.heappush(heap, (w, sv, to_v, eid))

        while heap:
            w, u, v, eid = heapq.heappop(heap)
            if v in visited:
                continue
            visited.add(v); remaining.discard(v)
            mst_eids.append(eid)
            total_w += w

            steps.append({
                'type': 'add_edge',
                'desc': f'เพิ่มเส้น ({u} — {v}), น้ำหนัก = {_fmt(w)}.  รวม = {_fmt(total_w)}',
                'visited': list(visited), 'mstEdges': list(mst_eids),
                'hlVerts': [u, v], 'hlEdges': [eid], 'totalWeight': total_w,
            })

            for to2, w2, eid2 in adj[v]:
                if to2 not in visited:
                    heapq.heappush(heap, (w2, v, to2, eid2))

    if start in v_ids:
        _prim_comp(start)
    while remaining:
        _prim_comp(min(remaining))

    label = 'MST' if len(roots) == 1 else 'Spanning Forest'
    steps.append({
        'type': 'result',
        'desc': (f'{label} เสร็จสมบูรณ์  น้ำหนักรวม = {_fmt(total_w)}'
                 + (f'  ({len(roots)} components)' if len(roots) > 1 else '')),
        'visited': list(visited), 'mstEdges': list(mst_eids),
        'hlVerts': list(visited), 'hlEdges': list(mst_eids),
        'totalWeight': total_w, 'componentCount': len(roots),
    })

    return {'steps': steps, 'totalWeight': total_w}


# ════════════════════════════════════════════════════════════
#  4. GRAPH COLORING  —  Welsh–Powell + Exact χ(G)
# ════════════════════════════════════════════════════════════

def _exact_chi(v_ids, adj):
    """Exact chromatic number via backtracking (|V| ≤ 15)."""
    n = len(v_ids)
    if n == 0:
        return 0
    cols = {}

    def bt(idx, k):
        if idx == n:
            return True
        v = v_ids[idx]
        used = {cols[nb] for nb in adj[v] if nb in cols}
        for c in range(k):
            if c not in used:
                cols[v] = c
                if bt(idx + 1, k):
                    return True
                del cols[v]
        return False

    for k in range(1, n + 1):
        cols.clear()
        if bt(0, k):
            return k
    return n


def welsh_powell(data):
    vertices = data['vertices']
    edges    = data['edges']
    if not vertices:
        return {'steps': [], 'colors': {}, 'chromaticNumber': 0}

    v_ids = [v['id'] for v in vertices]
    adj   = defaultdict(set)
    for e in edges:
        adj[e['from']].add(e['to'])
        adj[e['to']].add(e['from'])

    degs   = {v: len(adj[v]) for v in v_ids}
    sorted_v = sorted(v_ids, key=lambda x: degs[x], reverse=True)

    CNAMES = ['แดง','น้ำเงิน','เขียว','เหลือง','ม่วง','ส้ม',
              'ฟ้า','ชมพู','มะนาว','เขียวเข้ม','คราม','น้ำตาล',
              'บานเย็น','ทอง','เงิน']

    steps  = []
    colors = {}
    n_col  = 0

    steps.append({
        'type': 'init',
        'desc': ('เรียงจุดตาม Degree มาก→น้อย: '
                 + ', '.join(f'{v}(d={degs[v]})' for v in sorted_v)),
        'colors': {}, 'colorCount': 0,
        'hlVerts': sorted_v, 'hlEdges': [],
        'sortedVertices': sorted_v, 'degrees': degs,
    })

    for v in sorted_v:
        nb_cols = {colors[nb] for nb in adj[v] if nb in colors}
        c = 0
        while c in nb_cols:
            c += 1
        colors[v] = c
        if c >= n_col:
            n_col = c + 1
        cn = CNAMES[c] if c < len(CNAMES) else f'สี {c+1}'

        steps.append({
            'type': 'color',
            'desc': (f'ระบายจุด {v} (deg={degs[v]}) ด้วย {cn}'
                     f'  สีเพื่อนบ้าน: '
                     + (', '.join(CNAMES[x] if x < len(CNAMES) else f'สี {x+1}'
                                  for x in sorted(nb_cols)) or '—')),
            'colors': dict(colors), 'colorCount': n_col,
            'hlVerts': [v],
            'neighborVerts': [nb for nb in adj[v] if nb in colors],
            'hlEdges': [], 'currentVertex': v, 'assignedColor': c,
        })

    max_deg  = max(degs.values()) if degs else 0
    chi_exact = _exact_chi(v_ids, adj) if len(v_ids) <= 15 else None

    groups = defaultdict(list)
    for v, c in colors.items():
        cn = CNAMES[c] if c < len(CNAMES) else f'สี {c+1}'
        groups[cn].append(v)

    steps.append({
        'type': 'result',
        'desc': (f'ระบายสีเสร็จ: ใช้ {n_col} สี (Greedy).  Δ+1 = {max_deg+1}.'
                 + (f'  χ(G) = {chi_exact}' if chi_exact is not None else '')),
        'colors': dict(colors), 'colorCount': n_col,
        'hlVerts': v_ids, 'hlEdges': [],
        'colorGroups': dict(groups),
        'maxDegree': max_deg, 'chromaticNumber': chi_exact,
    })

    return {'steps': steps, 'colors': colors, 'colorCount': n_col,
            'chromaticNumber': chi_exact, 'colorGroups': dict(groups)}


# ════════════════════════════════════════════════════════════
#  5. NETWORK FLOW  —  Edmonds–Karp  (Ford–Fulkerson + BFS)
# ════════════════════════════════════════════════════════════

def edmonds_karp(data):
    vertices = data['vertices']
    edges    = data['edges']
    source   = data['source']
    sink     = data['sink']

    if not vertices or source == sink:
        return {'steps': [], 'maxFlow': 0}

    v_ids = [v['id'] for v in vertices]

    cap  = defaultdict(lambda: defaultdict(float))
    emap = defaultdict(list)
    for e in edges:
        cap[e['from']][e['to']] += e.get('weight', 1)
        emap[(e['from'], e['to'])].append(e.get('id', ''))

    res  = defaultdict(lambda: defaultdict(float))
    for u in cap:
        for v in cap[u]:
            res[u][v] = cap[u][v]

    flow = defaultdict(lambda: defaultdict(float))
    steps    = []
    max_flow = 0

    steps.append({
        'type': 'init',
        'desc': f'เริ่มต้น: source = {source}, sink = {sink}.  Flow ทั้งหมด = 0',
        'flow': {}, 'maxFlow': 0,
        'hlVerts': [source, sink], 'hlEdges': [], 'pathEdges': [],
        'augPath': [],
    })

    it = 0
    while True:
        it += 1
        parent = {source: None}
        vis = {source}
        q = deque([source])
        while q:
            u = q.popleft()
            if u == sink:
                break
            for v in sorted(set(list(res[u].keys()))):
                if v not in vis and res[u][v] > 1e-9:
                    vis.add(v)
                    parent[v] = u
                    q.append(v)

        if sink not in parent:
            break

        # Bottleneck
        bn = float('inf')
        c  = sink
        path = []
        while c is not None:
            path.append(c)
            p = parent[c]
            if p is not None:
                bn = min(bn, res[p][c])
            c = p
        path.reverse()

        pe = []
        for i in range(len(path)-1):
            eids = emap.get((path[i], path[i+1]), [])
            if eids:
                pe.append(eids[0])

        steps.append({
            'type': 'augmenting_path',
            'desc': (f'รอบ {it}: พบ Augmenting Path  '
                     f'{" → ".join(path)},  คอขวด = {_fmt(bn)}'),
            'flow': {}, 'maxFlow': max_flow,
            'hlVerts': path, 'hlEdges': pe, 'pathEdges': pe,
            'augPath': path, 'bottleneck': bn,
        })

        c = sink
        while c != source:
            p = parent[c]
            res[p][c] -= bn
            res[c][p] += bn
            flow[p][c] += bn
            flow[c][p] -= bn
            c = p
        max_flow += bn

        ef = {}
        for e in edges:
            eid = e.get('id', '')
            f   = flow[e['from']][e['to']]
            ca  = cap[e['from']][e['to']]
            ef[eid] = {'flow': max(0, f), 'capacity': ca}

        steps.append({
            'type': 'update_flow',
            'desc': f'อัปเดต Flow ตามเส้นทาง  คอขวด = {_fmt(bn)}.  Max Flow = {_fmt(max_flow)}',
            'flow': ef, 'maxFlow': max_flow,
            'hlVerts': [source, sink], 'hlEdges': pe, 'pathEdges': [],
            'augPath': path,
        })

    # Min-cut
    reach = set()
    q = deque([source]); reach.add(source)
    while q:
        u = q.popleft()
        for v in res[u]:
            if v not in reach and res[u][v] > 1e-9:
                reach.add(v)
                q.append(v)

    mc_eids = []
    for e in edges:
        if e['from'] in reach and e['to'] not in reach:
            mc_eids.append(e.get('id', ''))

    ef = {}
    for e in edges:
        eid = e.get('id', '')
        f   = flow[e['from']][e['to']]
        ca  = cap[e['from']][e['to']]
        ef[eid] = {'flow': max(0, f), 'capacity': ca}

    non_reach = sorted(set(v_ids) - reach)
    steps.append({
        'type': 'result',
        'desc': (f'Max Flow = {_fmt(max_flow)}.  '
                 f'Min Cut: {{{", ".join(sorted(reach))}}} | '
                 f'{{{", ".join(non_reach)}}}'),
        'flow': ef, 'maxFlow': max_flow,
        'hlVerts': list(reach), 'hlEdges': mc_eids,
        'pathEdges': [], 'minCutEdges': mc_eids,
        'reachable': sorted(reach), 'nonReachable': non_reach,
    })

    return {'steps': steps, 'maxFlow': max_flow,
            'minCutEdges': mc_eids, 'edgeFlows': ef}


# ════════════════════════════════════════════════════════════
#  HTTP SERVER
# ════════════════════════════════════════════════════════════

DISPATCH = {
    'analyze':  analyze_graph,
    'dijkstra': dijkstra,
    'prim':     prim,
    'coloring': welsh_powell,
    'flow':     edmonds_karp,
}


class Handler(http.server.BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        pass  # quiet

    def _cors(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'POST, GET, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')

    def do_OPTIONS(self):
        self.send_response(200)
        self._cors()
        self.end_headers()

    def do_GET(self):
        if self.path in ('/', '/index.html'):
            fpath = os.path.join(BASE_DIR, 'index.html')
            try:
                with open(fpath, 'r', encoding='utf-8') as f:
                    body = f.read().encode('utf-8')
                self.send_response(200)
                self.send_header('Content-Type', 'text/html; charset=utf-8')
                self.send_header('Content-Length', str(len(body)))
                self._cors()
                self.end_headers()
                self.wfile.write(body)
            except FileNotFoundError:
                self.send_response(404)
                self.end_headers()
                self.wfile.write(b'index.html not found')
        else:
            self.send_response(404)
            self.end_headers()

    def do_POST(self):
        length = int(self.headers.get('Content-Length', 0))
        raw    = self.rfile.read(length)
        data   = json.loads(raw.decode('utf-8'))
        action = data.get('action', '')

        try:
            fn     = DISPATCH.get(action)
            result = fn(data) if fn else {'error': f'Unknown action: {action}'}
        except Exception as exc:
            result = {'error': str(exc)}

        body = json.dumps(result, ensure_ascii=False).encode('utf-8')
        self.send_response(200)
        self.send_header('Content-Type', 'application/json; charset=utf-8')
        self.send_header('Content-Length', str(len(body)))
        self._cors()
        self.end_headers()
        self.wfile.write(body)


# ════════════════════════════════════════════════════════════
#  MAIN
# ════════════════════════════════════════════════════════════

def main():
    server = http.server.HTTPServer(('', PORT), Handler)
    import sys, io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    print()
    print('  +================================================+')
    print('  |  Graph Theory Lab -- Modern UI Edition          |')
    print('  |                                                 |')
    print(f'  |  http://localhost:{PORT}                         |')
    print('  |  Ctrl+C to stop                                 |')
    print('  +================================================+')
    print()

    threading.Timer(1.0, lambda: webbrowser.open(f'http://localhost:{PORT}')).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print('\n  🛑  Server stopped.')
        server.server_close()


if __name__ == '__main__':
    main()
