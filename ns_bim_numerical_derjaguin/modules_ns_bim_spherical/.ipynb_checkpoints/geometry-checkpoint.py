import numpy as np
import pandas as pd

################################################################################
# Functions for loading data

def get_node_links(file):
    # Assumes space in ' element_id' label and trailing comma in label list
    node_links = pd.read_csv(file)
    node_links.drop(labels=[' element_id'], axis=1, inplace=True)
    node_min = node_links.min().min()
    for i in node_links:
        for j in range(len(node_links[i])):
            node_links.at[j, i] -= (node_min - 1)
    node_links.dropna(axis=1, inplace=True)
    return node_links

def get_node_coords(file):
    # Assumes space in ' element_id' label and trailing comma in label list
    node_coords = pd.read_csv(file)
    node_coords.drop(labels=[' node_id'], axis=1, inplace=True)
    node_coords.columns = ['x_nm', 'y_nm', 'z_nm']
    return node_coords # Coordinates in [nm]

def get_potentials(file):
    # Assumes two rows (i.e. like two radii) in file
    potentials = pd.read_csv(file, header=None).transpose()
    potentials.columns = potentials.iloc[0]
    potentials = potentials[1:]
    potentials.columns = ['r1', 'r2']
    potentials.reset_index(inplace=True)
    potentials.drop(labels=['index','r2'], axis=1, inplace=True)

    # Convert potentials from mV to V
    for i in range(len(potentials['r1'])):
        potentials.at[i, 'r1'] *= 1.0e-3
    return potentials

################################################################################
# Functions for assisting in geometry calculations

def get_triangle_sides(vectors):
    """Get two sides of triangle whose corners are given as an array in vectors."""
    v1 = vectors[1] - vectors[0]
    v2 = vectors[2] - vectors[0]
    return np.array([v1, v2])

def get_triangle_area(node_vectors):
    """Compute area of triangle whose corners are given as an array in node_vectors."""
    sides = get_triangle_sides(node_vectors)
    v1 = sides[0]
    v2 = sides[1]
    return 0.5*np.linalg.norm(np.cross(v1, v2))

def get_protein_radius(elements_df):
    """Get protein radius."""
    return np.sqrt(elements_df.at[0, 'x1_nm']**2 + elements_df.at[0, 'y1_nm']**2 + \
                   elements_df.at[0, 'z1_nm']**2)

def get_characterstic_plane_point(elements_df, row, h):
    """Find the characteristic point of the plane separated a distance h
       from the element indicated by row in elements_df."""
    q_element = get_element_centroid_coords_vector(elements_df, row)
    n = get_element_normal_vector(elements_df, row)
    q_plane = q_element + h*n
    return q_plane

################################################################################
# Interface for retrieving boundary element data as numpy vectors/arrays

def get_element_node_coords_vectors(elements_df, row):
    """Retrieve array of vectors for the node coordinates of a single element."""
    xyz1 = [elements_df.at[row, 'x1_nm'], elements_df.at[row, 'y1_nm'],
            elements_df.at[row, 'z1_nm']]
    xyz2 = [elements_df.at[row, 'x2_nm'], elements_df.at[row, 'y2_nm'],
            elements_df.at[row, 'z2_nm']]
    xyz3 = [elements_df.at[row, 'x3_nm'], elements_df.at[row, 'y3_nm'],
            elements_df.at[row, 'z3_nm']]
    return np.array([xyz1, xyz2, xyz3])

def get_element_centroid_coords_vector(elements_df, row):
    """Retrieve vector of centroid coordinates for the element denoted by row."""
    xyz = [elements_df.at[row, 'x_cent_nm'], elements_df.at[row, 'y_cent_nm'],
           elements_df.at[row, 'z_cent_nm']]
    return np.array(xyz)

def get_element_normal_vector(elements_df, row):
    """Retrieve the normal vector for the element denoted by row."""
    xyz = [elements_df.at[row, 'n_x_nm'], elements_df.at[row, 'n_y_nm'],
           elements_df.at[row, 'n_z_nm']]
    return np.array(xyz)


################################################################################
# Functions for determining features of boundary elements

def set_element_potentials(node_links, potentials):
    for row in range(len(node_links['node_1_id'])):
        c1 = potentials.at[(node_links.at[row, 'node_1_id'] - 1), 'r1']
        c2 = potentials.at[(node_links.at[row, 'node_2_id'] - 1), 'r1']
        c3 = potentials.at[(node_links.at[row, 'node_3_id'] - 1), 'r1']
        node_links.at[row, 'potential_V'] = np.array([c1, c2, c3]).mean() # [mV]
    return node_links

def get_elements_node_coords_and_centroids(elements_df, node_coords_df):
    """Add node and centroid coordinates to the elements_df."""

    node_coords_array = node_coords_df.to_numpy()
    for row in range(len(elements_df)):
        node_1_id = elements_df.at[row, 'node_1_id']
        node_2_id = elements_df.at[row, 'node_2_id']
        node_3_id = elements_df.at[row, 'node_3_id']

        xyz1 = node_coords_array[node_1_id - 1]
        xyz2 = node_coords_array[node_2_id - 1]
        xyz3 = node_coords_array[node_3_id - 1]

        x_cent = np.mean(np.array([xyz1[0], xyz2[0], xyz3[0]]))
        y_cent = np.mean(np.array([xyz1[1], xyz2[1], xyz3[1]]))
        z_cent = np.mean(np.array([xyz1[2], xyz2[2], xyz3[2]]))

        elements_df.at[row, 'x1_nm'] = xyz1[0]
        elements_df.at[row, 'y1_nm'] = xyz1[1]
        elements_df.at[row, 'z1_nm'] = xyz1[2]
        elements_df.at[row, 'x2_nm'] = xyz2[0]
        elements_df.at[row, 'y2_nm'] = xyz2[1]
        elements_df.at[row, 'z2_nm'] = xyz2[2]
        elements_df.at[row, 'x3_nm'] = xyz3[0]
        elements_df.at[row, 'y3_nm'] = xyz3[1]
        elements_df.at[row, 'z3_nm'] = xyz3[2]
        elements_df.at[row, 'x_cent_nm'] = x_cent
        elements_df.at[row, 'y_cent_nm'] = y_cent
        elements_df.at[row, 'z_cent_nm'] = z_cent
    return

def get_elements_outward_unit_normal(elements_df):
    """Add outward unit normal vector coordinates to the elements_df."""

    for row in range(len(elements_df)):
        vectors = get_element_node_coords_vectors(elements_df, row)
        sides = get_triangle_sides(vectors)
        v1 = sides[0]
        v2 = sides[1]
        n = np.cross(v1, v2)
        n = n/np.linalg.norm(n)

        if np.sign(n[0]) == np.sign(elements_df.at[row, 'x_cent_nm']):
            pass
        else:
            n = -1.0*n

        elements_df.at[row, 'n_x_nm'] = n[0]
        elements_df.at[row, 'n_y_nm'] = n[1]
        elements_df.at[row, 'n_z_nm'] = n[2]
    return

def get_element_ids(elements_df):
    for row in range(len(elements_df)):
        elements_df.at[row, 'element_id'] = row
    return
################################################################################
# Master functions for generating the dictionary elements_dict, whose keys are
# the elements of protein_folders (eg. 'adh_pH_7'). These do not require
# optimization, as they are only performed once and take ~30 s.

def get_elem_and_node_dicts(nsbim_data_folder, protein_folders):
    elements_dict  = {}
    node_coords_dict = {}
    potentials_dict  = {}

    for i in protein_folders:
        data_dir = nsbim_data_folder + '/' + i
        node_links_file  = data_dir + '/rslt_node_linkages.csv'
        node_coords_file = data_dir + '/rslt_node_coordinates.csv'
        potentials_file  = data_dir + '/rslt_protein_potential_variable_radius.csv'

        elements_dict[i]  = get_node_links(node_links_file)
        node_coords_dict[i] = get_node_coords(node_coords_file)
        potentials_dict[i]  = get_potentials(potentials_file)

        elements_dict[i] = set_element_potentials(elements_dict[i], potentials_dict[i])

        # for j in range(len(elements_dict[i])):
        #     elements_dict[i].at[j, 'id'] = j
        # elements_dict[i] = elements_dict[i].astype({'id': 'int32'})

    return elements_dict, node_coords_dict

def fill_elem_dict(elements_dict, node_coords_dict, protein_folders):
    for i in range(len(protein_folders)):
        get_elements_node_coords_and_centroids(elements_dict[protein_folders[i]],
                                               node_coords_dict[protein_folders[i]])
        get_elements_outward_unit_normal(elements_dict[protein_folders[i]])
    return elements_dict


################################################################################
# Functions (v2) for hemisphere data, using the strategy of making a copy of
# elements_dict[key] for each CPU and applying vectorized functions.

def get_distance_helper(q_plane, n, x_cent, y_cent, z_cent):
    a_x = q_plane[0] - x_cent
    a_y = q_plane[1] - y_cent
    a_z = q_plane[2] - z_cent
    return np.abs(n[0]*a_x + n[1]*a_y + n[2]*a_z)

def get_distance(elements_df, row_plane):
    q_plane = get_characterstic_plane_point(elements_df, row_plane, 0.0)
    n       = get_element_normal_vector(elements_df, row_plane)

    elements_df['0_dist_nm'] = get_distance_helper(q_plane, n,
    elements_df['x_cent_nm'].values, elements_df['y_cent_nm'].values,
    elements_df['z_cent_nm'].values)
    return

def get_area_helper(n, q_plane, x1, y1, z1, x2, y2, z2, x3, y3, z3):
    a = n[0]
    b = n[1]
    c = n[2]
    d = q_plane[0]
    e = q_plane[1]
    f = q_plane[2]

    # Projected points px, py, pz
    t = (a*(d-x1) + b*(e-y1) + c*(f-z1))/(a**2 + b**2 + c**2)
    (px1, py1, pz1) = ((x1 + t*a), (y1 + t*b), (z1 + t*c))

    t = (a*(d-x2) + b*(e-y2) + c*(f-z2))/(a**2 + b**2 + c**2)
    (px2, py2, pz2) = ((x2 + t*a), (y2 + t*b), (z2 + t*c))

    t = (a*(d-x3) + b*(e-y3) + c*(f-z3))/(a**2 + b**2 + c**2)
    (px3, py3, pz3) = ((x3 + t*a), (y3 + t*b), (z3 + t*c))

    # Triangle side vectors
    (sx1, sy1, sz1) = (px2 - px1, py2 - py1, pz2 - pz1)
    (sx2, sy2, sz2) = (px3 - px1, py3 - py1, pz3 - pz1)

    # Cross product of the triangle sides
    cx = sy1*sz2 - sz1*sy2
    cy = sz1*sx2 - sx1*sz2
    cz = sx1*sy2 - sy1*sx2
    norm = np.sqrt(cx**2 + cy**2 + cz**2)
    return 0.5*norm

def get_area(elements_df, row_plane):
    n       = get_element_normal_vector(elements_df, row_plane)
    q_plane = get_characterstic_plane_point(elements_df, row_plane, 0.0)

    elements_df['area_nm2'] = get_area_helper(n, q_plane,
    elements_df['x1_nm'].values, elements_df['y1_nm'].values, elements_df['z1_nm'].values,
    elements_df['x2_nm'].values, elements_df['y2_nm'].values, elements_df['z2_nm'].values,
    elements_df['x3_nm'].values, elements_df['y3_nm'].values, elements_df['z3_nm'].values)
    return

def cut_list_helper(cut_off, distance):
    return np.greater(distance, cut_off) # cut where True

def cut(elements_df):
    cut_off  = get_protein_radius(elements_df)
    cut_list = cut_list_helper(cut_off, elements_df['0_dist_nm'])

    cut_indeces = []
    for i in range(len(cut_list)):
        if cut_list[i] == True:
            cut_indeces.append(i)

    elements_df.drop(labels=cut_indeces, axis='index', inplace=True)
    return

def get_orientation_df(elements_df, row_plane):
    orientation_df = elements_df.copy()
    get_distance(orientation_df, row_plane)
    cut(orientation_df)
    get_area(orientation_df, row_plane)
    return orientation_df
























################################################################################
# Boneyard

def get_element_areas(elements_df, df_hemi):
    """Get element areas (no projection into plane)."""
    areas  = []
    for i in range(len(df_hemi)):
        row_element = int(df_hemi.at[i, 'element_row_num'])
        nodes = get_element_node_coords_vectors(elements_df, row_element)
        areas.append(get_triangle_area(nodes))  # [nm^2]
    df_hemi['area_nm2'] = areas
    return df_hemi

# Functions (v1) for generating df_hemi. These are suboptimal and not used.
def get_distance_to_plane(elements_df, row_plane, row_centroid, h):
    """Compute distance from centroid to plane defined by row_plane and h.
       Perhaps make this part of the main program to remove repeat calculations
       of q_plane and n."""
    q_plane    = get_characterstic_plane_point(elements_df, row_plane, h)
    n          = get_element_normal_vector(elements_df, row_plane)
    p_centroid = get_element_centroid_coords_vector(elements_df, row_centroid)
    distance   = np.dot((q_plane - p_centroid), n)
    return distance

def get_nodes_plane_projection(elements_df, row_plane, row_nodes, h):
    """Get projections of nodes for the element row_nodes into the plane defined
       by row_plane and h."""
    n            = get_element_normal_vector(elements_df, row_plane)
    q_plane      = get_characterstic_plane_point(elements_df, row_plane, h)
    node_vectors = get_element_node_coords_vectors(elements_df, row_nodes)

    a = n[0]
    b = n[1]
    c = n[2]
    d = q_plane[0]
    e = q_plane[1]
    f = q_plane[2]

    node_projections = []
    for node in node_vectors:
        x = node[0]
        y = node[1]
        z = node[2]

        t = (a*(d-x) + b*(e-y) + c*(f-z))/(a**2 + b**2 + c**2)
        projection = [(x + t*a), (y + t*b), (z + t*c)]
        node_projections.append(projection)
    return np.array(node_projections)

def get_hemisphere_elements(elements_df, row_plane):
    """Identify the elements in the hemisphere facing the plane defined by
       row_plane in elements_df. Set up a DataFrame to hold their data."""
    radius = get_protein_radius(elements_df)
    df_hemi = pd.DataFrame(columns=['element_row_num'])

    element_row_nums = []
    for row in range(len(elements_df)):
        d = get_distance_to_plane(elements_df, row_plane, row_centroid=row, h=0.0)
        if d < radius:
            element_row_nums.append(row)
    df_hemi = pd.DataFrame(element_row_nums, columns=['element_row_num'])
    return df_hemi

def get_hemisphere_potentials(elements_df, df_hemi):
    """Get potentials for hemisphere elements."""
    df_hemi['potential_V'] = np.nan
    for i in range(len(df_hemi)):
        row_element = int(df_hemi.at[i, 'element_row_num'])
        df_hemi.at[i, 'potential_V'] = elements_df.at[row_element, 'potential_V']
    return df_hemi

def get_hemisphere_distances(elements_df, row_plane, h, df_hemi):
    """Get the distances between the hemisphere centroid points and plane defined
       by row_plane in elements_df and h."""
    distances = []
    for i in range(len(df_hemi)):
        row_element = int(df_hemi.at[i, 'element_row_num'])
        distances.append(get_distance_to_plane(elements_df, row_plane, row_element, h))
    df_hemi['distance_nm'] = distances
    return df_hemi

def get_projected_areas(elements_df, row_plane, df_hemi):
    """Get element projected areas into the plane defined by row_plane."""
    areas = []
    for i in range(len(df_hemi)):
        row_element = int(df_hemi.at[i, 'element_row_num'])
        projected_nodes = get_nodes_plane_projection(elements_df, row_plane,
                                                     row_element, h=0.0)
        areas.append(get_triangle_area(projected_nodes))    # [nm^2]
    df_hemi['proj_area_nm2'] = areas
    return df_hemi

def get_df_hemi(elements_df, surf_element, z_0):
    df_hemi = get_hemisphere_elements(elements_df, surf_element)
    df_hemi = get_hemisphere_potentials(elements_df, df_hemi)
    df_hemi = get_hemisphere_distances(elements_df, surf_element, z_0, df_hemi)
    df_hemi = get_projected_areas(elements_df, surf_element, df_hemi)
    return df_hemi
