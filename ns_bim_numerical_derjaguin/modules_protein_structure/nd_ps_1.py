
################################################################################
# Get protein surface data

def get_apbs_data(dir_apbs, protein_name, file):
    df = pd.read_csv(dir_apbs + protein_name + '/surf_pot_' + file + '.csv',
                    names=['x_A', 'y_A', 'z_A', 'phi_V'])
    df['phi_V'] *= constants().kT/constants().e0
    return df

def get_tabi_data(dir_tabi, protein_name):
    df = pd.read_csv(dir_tabi + protein_name + '/rslt_surface_potential.csv',
                    names=['node_index', 'x_A', 'y_A', 'z_A', 'norm_x', 'norm_y',
                            'norm_z', 'phi_V', 'norm_phi'])
    df['phi_V'] *= 0.0433699
    return df

def get_link(dir_link, protein_name):
    df = pd.read_csv(dir_link + protein_name + '/rslt_node_links.csv', names=['node_1',
                    'node_2', 'node_3'])
    for col in df.columns:
        df[col] -= 1
    return df

def get_elements(df_link, df_data):
    for ind, contents in df_link.iterrows():
        df_link.at[ind, 'x1_A'] = df_data.at[contents['node_1'], 'x_A']
        df_link.at[ind, 'y1_A'] = df_data.at[contents['node_1'], 'y_A']
        df_link.at[ind, 'z1_A'] = df_data.at[contents['node_1'], 'z_A']
        df_link.at[ind, 'x2_A'] = df_data.at[contents['node_2'], 'x_A']
        df_link.at[ind, 'y2_A'] = df_data.at[contents['node_2'], 'y_A']
        df_link.at[ind, 'z2_A'] = df_data.at[contents['node_2'], 'z_A']
        df_link.at[ind, 'x3_A'] = df_data.at[contents['node_3'], 'x_A']
        df_link.at[ind, 'y3_A'] = df_data.at[contents['node_3'], 'y_A']
        df_link.at[ind, 'z3_A'] = df_data.at[contents['node_3'], 'z_A']

        phi_1 = df_data.at[contents['node_1'], 'phi_V']
        phi_2 = df_data.at[contents['node_2'], 'phi_V']
        phi_3 = df_data.at[contents['node_3'], 'phi_V']
        df_link.at[ind, 'phi_V'] = np.array([phi_1, phi_2, phi_3]).mean()

        df_link.at[ind, 'x_cent_A'] = np.mean(np.array([df_link.at[ind, 'x1_A'],
                                                        df_link.at[ind, 'x2_A'],
                                                        df_link.at[ind, 'x3_A']]))
        df_link.at[ind, 'y_cent_A'] = np.mean(np.array([df_link.at[ind, 'y1_A'],
                                                        df_link.at[ind, 'y2_A'],
                                                        df_link.at[ind, 'y3_A']]))
        df_link.at[ind, 'z_cent_A'] = np.mean(np.array([df_link.at[ind, 'z1_A'],
                                                        df_link.at[ind, 'z2_A'],
                                                        df_link.at[ind, 'z3_A']]))
    return df_link

def translate_prot_centroid(df_prot):
    n_elements = len(df_prot['x_cent_A'])

    prot_x_cent = sum(df_prot['x_cent_A'])/n_elements
    prot_y_cent = sum(df_prot['y_cent_A'])/n_elements
    prot_z_cent = sum(df_prot['z_cent_A'])/n_elements

    df_prot['x_cent_A'] -= prot_x_cent
    df_prot['y_cent_A'] -= prot_y_cent
    df_prot['z_cent_A'] -= prot_z_cent
    return

def get_df_elem_apbs(dir_link, dir_apbs, protein_name, file):
    df_link = get_link(dir_link, protein_name)
    df_data = get_apbs_data(dir_apbs, protein_name, file)
    df_prot = get_elements(df_link, df_data)
    translate_prot_centroid(df_prot)
    label_edges(df_prot)
    return df_prot

def get_df_elem_tabi(dir_link, dir_apbs, protein_name):
    df_link = get_link(dir_link, protein_name)
    df_data = get_tabi_data(dir_tabi, protein_name)
    df_prot = get_elements(df_link, df_data)
    translate_prot_centroid(df_prot)
    label_edges(df_prot)
    return df_prot

def label_edges(df):
    df['edge_1'] = np.nan
    df['edge_2'] = np.nan
    df['edge_3'] = np.nan

    df['edge_1'] = df['edge_1'].astype('object')
    df['edge_2'] = df['edge_2'].astype('object')
    df['edge_3'] = df['edge_3'].astype('object')

    for ind, content in df.iterrows():
        e1 = [content['node_1'], content['node_2']]
        e2 = [content['node_1'], content['node_3']]
        e3 = [content['node_2'], content['node_3']]
        e1.sort()
        e2.sort()
        e3.sort()

        df.at[ind, 'edge_1'] = tuple(e1)
        df.at[ind, 'edge_2'] = tuple(e2)
        df.at[ind, 'edge_3'] = tuple(e3)
    return


################################################################################
# Helper functions for performing 3D geometry calculations

def get_unit_vector(theta, phi):
    """Find unit vector pointed in the theta, phi direction."""
    a = np.sin(theta)*np.cos(phi)
    b = np.sin(theta)*np.sin(phi)
    c = np.cos(theta)
    return [a, b, c]

def get_distance_helper(n, q_plane, x, y, z):
    """Find the distance between a point (x, y, z) and a plane (defined by
    the unit normal n and the point q_plane)."""
    a  = n[0]
    b  = n[1]
    c  = n[2]
    x0 = q_plane[0]
    y0 = q_plane[1]
    z0 = q_plane[2]
    d  = a*x0 + b*y0 + c*z0
    distance = a*x + b*y + c*z - d
    return distance

def get_projection_helper(n, q_plane, x, y, z):
    """Find the boundary element's projection into the resin plane."""
    a = n[0]
    b = n[1]
    c = n[2]
    d = q_plane[0]
    e = q_plane[1]
    f = q_plane[2]

    # Projected points px, py, pz
    t = (a*(d-x) + b*(e-y) + c*(f-z))/(a**2 + b**2 + c**2)
    (px, py, pz) = ((x + t*a), (y + t*b), (z + t*c))
    return (px, py, pz)

def translate_to_origin_helper(n, x, y, z):
    (nx, ny, nz) = n
    q_plane  = [x, y, z]
    distance = get_distance_helper(n, q_plane, 0, 0, 0)
    (tx, ty, tz) = (x + distance*nx, y + distance*ny, z + distance*nz)
    return (tx, ty, tz)

def rotate_to_xy_plane_helper(n, x, y, z):
    """Rotate point (x, y, z) on plane with unit normal n into the xy plane."""
    a = n[0]
    b = n[1]
    c = n[2]

    r11 = c + (b**2)*(1.0-c)/(a**2 + b**2)
    r12 = -1.0*a*b*(1.0-c)/(a**2 + b**2)
    r13 = -1.0*a
    r21 = r12
    r22 = c + (a**2)*(1.0-c)/(a**2 + b**2)
    r23 = -1.0*b
    r31 = a
    r32 = b
    r33 = c

    tx = r11*x + r12*y + r13*z
    ty = r21*x + r22*y + r23*z
    tz = r31*x + r32*y + r33*z
    return (tx, ty, tz)

def get_side_distance(x1, y1, z1, x2, y2, z2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)

def get_triangle_area_helper(px1, py1, pz1, px2, py2, pz2, px3, py3, pz3):
    # Triangle side vectors
    (sx1, sy1, sz1) = (px2 - px1, py2 - py1, pz2 - pz1)
    (sx2, sy2, sz2) = (px3 - px1, py3 - py1, pz3 - pz1)

    # Cross product of the triangle sides
    cx = sy1*sz2 - sz1*sy2
    cy = sz1*sx2 - sx1*sz2
    cz = sx1*sy2 - sy1*sx2
    norm = np.sqrt(cx**2 + cy**2 + cz**2)
    return 0.5*norm

# def store_xy_vertices_helper(x1, y1, x2, y2, x3, y3):
#     return [[x1, y1], [x2, y2], [x3, y3]]

################################################################################
# 3D geometry calculations on full protein structure

def get_distance(n, q_plane, protein):
    protein.df['distance_A'] = get_distance_helper(n, q_plane,
    protein.df['x_cent_A'].values, protein.df['y_cent_A'].values,
    protein.df['z_cent_A'].values)
    return

# # Currently unused
# def get_resin_projection_coords(n, q_plane, df):
#     (df['px_A'], df['py_A'], df['pz_A']) =\
#     get_projection_helper(n, q_plane, df['x_A'].values, df['y_A'].values,
#                           df['z_A'].values)
#     return
#
# def translate_projection_to_origin_coords(n, df):
#     (df['px_A'], df['py_A'], df['pz_A']) =\
#     translate_to_origin_helper(n, df['px_A'].values, df['py_A'].values,
#                                df['pz_A'].values)
#     return
#
# def rotate_to_xy_plane_coords(n, df):
#     (df['px_A'], df['py_A'], df['pz_A']) =\
#     rotate_to_xy_plane_helper(n, df['px_A'].values, df['py_A'].values,
#                               df['pz_A'].values)
#     return

# Could reduce code reduncancy by using the functions above, with modification
# for column names - would require passing columns directly
def get_resin_projection(n, q_plane, protein):
    (protein.df['px1_A'], protein.df['py1_A'], protein.df['pz1_A']) =\
    get_projection_helper(n, q_plane, protein.df['x1_A'].values,
                          protein.df['y1_A'].values, protein.df['z1_A'].values)

    (protein.df['px2_A'], protein.df['py2_A'], protein.df['pz2_A']) =\
    get_projection_helper(n, q_plane, protein.df['x2_A'].values,
                          protein.df['y2_A'].values, protein.df['z2_A'].values)

    (protein.df['px3_A'], protein.df['py3_A'], protein.df['pz3_A']) =\
    get_projection_helper(n, q_plane, protein.df['x3_A'].values,
                          protein.df['y3_A'].values, protein.df['z3_A'].values)
    return

def translate_projection_to_origin(n, protein):
    (protein.df['px1_A'], protein.df['py1_A'], protein.df['pz1_A']) =\
    translate_to_origin_helper(n, protein.df['px1_A'].values,
                          protein.df['py1_A'].values, protein.df['pz1_A'].values)

    (protein.df['px2_A'], protein.df['py2_A'], protein.df['pz2_A']) =\
    translate_to_origin_helper(n, protein.df['px2_A'].values,
                          protein.df['py2_A'].values, protein.df['pz2_A'].values)

    (protein.df['px3_A'], protein.df['py3_A'], protein.df['pz3_A']) =\
    translate_to_origin_helper(n, protein.df['px3_A'].values,
                          protein.df['py3_A'].values, protein.df['pz3_A'].values)
    return

def rotate_to_xy_plane(n, protein):
    (protein.df['px1_A'], protein.df['py1_A'], protein.df['pz1_A']) =\
    rotate_to_xy_plane_helper(n, protein.df['px1_A'].values,
                          protein.df['py1_A'].values, protein.df['pz1_A'].values)

    (protein.df['px2_A'], protein.df['py2_A'], protein.df['pz2_A']) =\
    rotate_to_xy_plane_helper(n, protein.df['px2_A'].values,
                          protein.df['py2_A'].values, protein.df['pz2_A'].values)

    (protein.df['px3_A'], protein.df['py3_A'], protein.df['pz3_A']) =\
    rotate_to_xy_plane_helper(n, protein.df['px3_A'].values,
                          protein.df['py3_A'].values, protein.df['pz3_A'].values)
    return

def get_p_side_lengths(protein):
    protein.df['p_side_1_A'] = get_side_distance(protein.df['px1_A'].values,
    protein.df['py1_A'].values, protein.df['pz1_A'].values, protein.df['px2_A'].values,
    protein.df['py2_A'].values, protein.df['pz2_A'].values)

    protein.df['p_side_2_A'] = get_side_distance(protein.df['px1_A'].values,
    protein.df['py1_A'].values, protein.df['pz1_A'].values, protein.df['px3_A'].values,
    protein.df['py3_A'].values, protein.df['pz3_A'].values)

    protein.df['p_side_3_A'] = get_side_distance(protein.df['px3_A'].values,
    protein.df['py3_A'].values, protein.df['pz3_A'].values, protein.df['px2_A'].values,
    protein.df['py2_A'].values, protein.df['pz2_A'].values)
    return

def store_xy_vertices(protein):
    protein.df['vertices_A'] = np.nan
    protein.df['vertices_A'] = protein.df['vertices_A'].astype('object')
    protein.df['poly'] = np.nan
    protein.df['poly'] = protein.df['poly'].astype('object')

    for ind, content in protein.df.iterrows():
        protein.df.at[ind, 'vertices_A'] = [[content['px1_A'],
        content['py1_A']], [content['px2_A'], content['py2_A']], [content['px3_A'],
        content['py3_A']]]
        protein.df.at[ind, 'poly'] = geometry.Polygon(protein.df.at[ind, 'vertices_A'])
    return

# def store_xy_vertices_2(protein):
#     protein.df['vertices_A'] = np.nan
#     protein.df['vertices_A'] = protein.df['vertices_A'].astype('object')
#     p_vert_zip = zip()
#     protein.df['vertices_A'] = store_xy_vertices_helper(protein.df['px1_A'],
#         protein.df['py1_A'], protein.df['px2_A'], protein.df['py2_A'],
#         protein.df['px3_A'], protein.df['py3_A'])
#     return

################################################################################
# Functions for 2D geometry, etc.
def get_max_p_side_length(protein):
    a = max(protein.df['p_side_1_A'])
    b = max(protein.df['p_side_2_A'])
    c = max(protein.df['p_side_3_A'])
    return max(a, b, c)

def get_areas(protein):
    max_side = get_max_p_side_length(protein)
    i_cnt = 0

    for ind, contents in protein.df.iterrows():
        subj_dist = contents['distance_A']
        subj_poly = geometry.Polygon(contents['vertices_A'])

        # Find the box bounding the subject triangle
        x_min = min(contents['px1_A'], contents['px2_A'], contents['px3_A'])
        x_max = max(contents['px1_A'], contents['px2_A'], contents['px3_A'])
        y_min = min(contents['py1_A'], contents['py2_A'], contents['py3_A'])
        y_max = max(contents['py1_A'], contents['py2_A'], contents['py3_A'])

        # Extend the bounding box to catch all crossing triangles
        x_adj = 1.001*(max_side - (x_max - x_min))
        y_adj = 1.001*(max_side - (y_max - y_min))
        x_min -= x_adj
        x_max += x_adj
        y_min -= y_adj
        y_max += y_adj

        # Find the triangles which (1) have at least one point in the extended bounding box and (2) are closer to the resin than the subject
        a = protein.df[(((x_min < protein.df['px1_A']) & (protein.df['px1_A'] < x_max) & (y_min < protein.df['py1_A']) & (protein.df['py1_A'] < y_max)) |\
                        ((x_min < protein.df['px2_A']) & (protein.df['px2_A'] < x_max) & (y_min < protein.df['py2_A']) & (protein.df['py2_A'] < y_max)) |\
                        ((x_min < protein.df['px3_A']) & (protein.df['px3_A'] < x_max) & (y_min < protein.df['py3_A']) & (protein.df['py3_A'] < y_max))) &\
                         (protein.df['distance_A'] < subj_dist)]

        m = []
        for a_ind in a.index:
            m.append(geometry.Polygon(a.at[a_ind, 'vertices_A']))

        interference = ops.unary_union(m)
        try:
            protein.df.at[ind, 'area_A2'] = subj_poly.difference(interference).area
        except:
            i_cnt += 1
            print(i_cnt, 'Issue', ind)
            protein.df.at[ind, 'area_A2'] = 0.0

        if ind%1000 == 0:
            print(ind)
    return



################################################################################
# Driver functions

def main_projection(n, q_plane, protein):
    get_distance(n, q_plane, protein)

    get_resin_projection(n, q_plane, protein)
    translate_projection_to_origin(n, protein)
    rotate_to_xy_plane(n, protein)
    # round_p_coords(protein)
    get_p_side_lengths(protein)

    store_xy_vertices(protein)
    # store_xy_vertices_2(protein)
    return

def get_min_helper(x1, x2, x3):
    a = np.minimum(x1, x2)
    return np.minimum(a, x3)

def get_max_helper(x1, x2, x3):
    a = np.maximum(x1, x2)
    return np.maximum(a, x3)

def get_obstructions_helper(x_min, x_max, y_min, y_max, subj_dist, protein):
    a = protein.df[(((x_min < protein.df['px1_A']) & (protein.df['px1_A'] < x_max) & (y_min < protein.df['py1_A']) & (protein.df['py1_A'] < y_max)) |\
                    ((x_min < protein.df['px2_A']) & (protein.df['px2_A'] < x_max) & (y_min < protein.df['py2_A']) & (protein.df['py2_A'] < y_max)) |\
                    ((x_min < protein.df['px3_A']) & (protein.df['px3_A'] < x_max) & (y_min < protein.df['py3_A']) & (protein.df['py3_A'] < y_max))) &\
                     (protein.df['distance_A'] < subj_dist)]
    return a

def get_interference_helper(protein):
    protein.df['interference'] = np.nan
    protein.df['interference'] = protein.df['interference'].astype('object')

    for ind, content in protein.df.iterrows():
        protein.df.at[ind, 'interference'] = ops.unary_union(content['obstructions']['poly'])
    return

def get_areas_2(protein):
    max_side = get_max_p_side_length(protein)

    protein.df['x_min'] = get_min_helper(protein.df['px1_A'], protein.df['px2_A'], protein.df['px3_A'])
    protein.df['x_max'] = get_max_helper(protein.df['px1_A'], protein.df['px2_A'], protein.df['px3_A'])
    protein.df['y_min'] = get_min_helper(protein.df['py1_A'], protein.df['py2_A'], protein.df['py3_A'])
    protein.df['y_max'] = get_max_helper(protein.df['py1_A'], protein.df['py2_A'], protein.df['py3_A'])

    # Extend the bounding box to catch all crossing triangles
    x_adj = 1.001*(max_side - (protein.df['x_max'] - protein.df['x_min']))
    y_adj = 1.001*(max_side - (protein.df['y_max'] - protein.df['y_min']))
    protein.df['x_min'] -= x_adj
    protein.df['x_max'] += x_adj
    protein.df['y_min'] -= y_adj
    protein.df['y_max'] += y_adj

    protein.df['obstructions'] = np.nan
    protein.df['obstructions'] = protein.df['poly'].astype('object')
    temp = np.vectorize(get_obstructions_helper)
    protein.df['obstructions'] = temp(protein.df['x_min'], protein.df['x_max'],
        protein.df['y_min'], protein.df['y_max'], protein.df['distance_A'], protein)

    # m = []
    # for a_ind in a.index:
    #     m.append(geometry.Polygon(a.at[a_ind, 'vertices_A']))
    #
    # interference = ops.unary_union(m)
    # try:
    #     protein.df.at[ind, 'area_A2'] = subj_poly.difference(interference).area
    # except:
    #     i_cnt += 1
    #     print(i_cnt, 'Issue', ind)
    #     protein.df.at[ind, 'area_A2'] = 0.0
    #
    # if ind%1000 == 0:
    #     print(ind)
    return




################################################################################




def store_polygons(protein):
    protein.df['poly'] = np.nan
    protein.df['poly'] = protein.df['poly'].astype('object')

    for ind, content in protein.df.iterrows():
        protein.df.at[ind, 'poly'] = geometry.Polygon(content['vertices_A'])
    return




































################################################################################
# Boneyard


# def get_df_elem_apbs(dir_link, dir_apbs, protein_name, file):
#     df_link = get_link(dir_link, protein_name)
#     df_data = get_apbs_data(dir_apbs, protein_name, file)
#     df_prot = get_elements(df_link, df_data)
#     df_data.drop(labels=['phi_V'], axis='columns', inplace=True)
#     return (df_prot, df_data)
#
# def get_df_elem_tabi(dir_link, dir_apbs, protein_name):
#     df_link = get_link(dir_link, protein_name)
#     df_data = get_tabi_data(dir_tabi, protein_name)
#     df_prot = get_elements(df_link, df_data)
#     df_data.drop(labels=['node_index', 'norm_x', 'norm_y', 'norm_z', 'phi_V',
#                          'norm_phi'], axis='columns', inplace=True)
#     return (df_prot, df_data)



# def refine_vertices(protein):
#     clip_verts = protein.df.at[0, 'vertices_A']
#     protein.df['vertices_A'] = clip(protein.df['vertices_A'], clip_verts)
#     return

# # Not working
# def get_epicenter(x1, x2, x3, y1, y2, y3):
#     y = ((x1-x3)*(y2*2 - y1**2 + x2**2 - x1**2) + (x2-x1)*(x3**2 - x1**2 + y3**2 - y1**2))/\
#         (2.0*((y3-y1)*(x2-x1) - (x1-x3)*(y1-y2)))
#     x = y*(y1-y2)/(x2-x1) + (y2*2 - y1**2 + x2**2 - x1**2)/(2.0*(x2-x1))
#     return (x, y)

# def sort_vertices_counterclockwise(array):
#     s = np.where(array != 0)
#     yr = s[0].astype(np.uint8)
#     xc = s[1].astype(np.uint8)
#     center_xc = np.sum(xc)/xc.shape
#     center_yr = np.sum(yr)/yr.shape
#     theta = np.arctan2(yr-center_yr, xc-center_xc) * 180 / np.pi
#     indices = np.argsort(theta)
#     x = xc[indices]
#     y = yr[indices]
#     return x, y

# ################################################################################
# # 2D geometry calculations in the xy plane
#
# def list_unique_edges(seq):
#     seen = set()
#     seen_add = seen.add
#     # adds all elements it doesn't know yet to seen and all other to seen_twice
#     seen_twice = set( x for x in seq if x in seen or seen_add(x) )
#     seen_once  = set( x for x in seq if x not in seen_twice)
#     return list(seen_once)
#
# def find_unique_polygons(edges):
#     # Get unique node IDs
#     k = [inner for outer in edges for inner in outer]
#     nodes = []
#     [nodes.append(x) for x in k if x not in nodes]
#
#     # Make dictionary where key is the node ID, values are the indeces
#     # of edge elements which contain the node
#     m = {}
#     for x in nodes:
#         m[x] = [edges.index(outer) for outer in edges if x in outer]
#
#     # Trace out the node connections of each polygon
#     # z is the number of nodes accounted for; eind is the edge element index
#     poly_nodes = []
#     z = 0
#     eind = 0
#
#     while z < len(m):
#         res_eind = []
#         res_node = []
#
#         node = edges[eind][0]
#         k = [inner for outer in poly_nodes for inner in outer]
#
#         # Starting with an unused edge
#         while node in k:
#             eind += 1
#             node = edges[eind][0]
#
#         res_eind.append(eind)
#         res_node.append(node)
#
#         # Trace the node connections for the present polygon
#         cnt = 0
#         while cnt < len(m):
#             for x in m[node]:
#                 if x not in res_eind:
#                     eind = x
#
#             for x in edges[eind]:
#                 if x not in res_node:
#                     node = x
#
#             if node in res_node:
#                 break
#             else:
#                 res_eind.append(eind)
#                 res_node.append(node)
#                 cnt += 1
#
#         # Store the present polygon's node connectivity
#         poly_nodes.append(res_node)
#         for x in poly_nodes:
#             z += len(x)
#
#     return poly_nodes
#
# def find_polygon_vertices(poly_nodes, df_coords):
#     poly_vert = []
#     for x in poly_nodes:
#         res = []
#         for node in x:
#             res.append([df_coords.at[node, 'px_A'], df_coords.at[node, 'py_A']])
#         poly_vert.append(res)
#     return poly_vert
#
# def encompass_test(subject_vert, poly_vert):
#     poly_in = Polygon(subject_vert)
#     encompass_test = []
#     for x in poly_vert:
#         poly_out = Polygon(x)
#         encompass_test.append(poly_out.contains(poly_in))
#     if True in encompass_test:
#         return True
#     else:
#         return False
#
# def clip(subjectPolygon, clipPolygon):
#     # Adapted from the source code at:
#     # https://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Python
#     def inside(p):
#         return(cp2[0]-cp1[0])*(p[1]-cp1[1]) > (cp2[1]-cp1[1])*(p[0]-cp1[0])
#
#     def computeIntersection():
#         dc = [ cp1[0] - cp2[0], cp1[1] - cp2[1] ]
#         dp = [ s[0] - e[0], s[1] - e[1] ]
#         n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
#         n2 = s[0] * e[1] - s[1] * e[0]
#         n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0])
#         return [(n1*dp[0] - n2*dc[0]) * n3, (n1*dp[1] - n2*dc[1]) * n3]
#
#     outputList = subjectPolygon
#     cp1 = clipPolygon[-1]
#
#     cnt = 0
#     for clipVertex in clipPolygon:
#         cp2 = clipVertex
#         inputList = outputList
#         outputList = []
#         s = inputList[-1]
#         # print(inputList)
#
#         for subjectVertex in inputList:
#             e = subjectVertex
#             if inside(e):
#                 if not inside(s):
#                     outputList.append(computeIntersection())
#                 outputList.append(e)
#             elif inside(s):
#                 outputList.append(computeIntersection())
#             s = e
#         cp1 = cp2
#         print(outputList)
#
#     # Modification to remove duplicates from vertices list
#     res = []
#     [res.append(x) for x in outputList if x not in res]
#     return res
#
# ################################################################################

# def round_p_coords(protein):
#     keys = ['px1_A', 'py1_A', 'pz1_A', 'px2_A', 'py2_A', 'pz2_A', 'px3_A',
#             'py3_A', 'pz3_A']
#     decimals = {}
#     for key in keys:
#         decimals[key] = 4
#     protein.df.round(decimals)
#     return
