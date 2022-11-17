def read_file(file_path):
    # read the file
    with open(file_path) as f:
        # Get all file lines
        lines = f.readlines()
        timestep_list = []
        for i in range(0, len(lines)):
            textline = lines[i].strip().split()
            if "TIMESTEP" in textline:
                #read lattice parameters in case of rearrange by pbc
                x_lattice = get_cell(lines[i + 5])[1]
                y_lattice = get_cell(lines[i + 6])[1]
                z_lattice = get_cell(lines[i + 7])[1] - get_cell(lines[i + 7])[0]
                latpar = [float(x_lattice), float(y_lattice), float(z_lattice)]
                # reboot temporal list
                if i > 100:
                    atoms_list = [pb_list, i_list, latpar]
                    timestep_list.append(atoms_list)
                pb_list = []
                i_list = []
            elif len(textline) > 1 and textline[1] == '1':
                # read line
                pb_line = lines[i]
                # get coords
                pb_coord = get_coord(pb_line)
                # use boundary conditions to replicate atoms in the lateral faces of the supercell
                pb_bound_cond = replica(pb_coord, latpar)
                # merge identical atoms  (by periodic boundary conditions) in the same list
                pb_identical_atoms = [pb_bound_cond]
                pb_list.extend(pb_identical_atoms)
            elif len(textline) > 1 and textline[1] == '2':
                i_line = lines[i]
                i_coord = get_coord(i_line)
                i_list.append(i_coord)

    return timestep_list


def replica(pb_atom, latpar):
    # make replicas of pb atoms in xz and yz planes due to periodic boundary conditions
    replicas = []
    #save old coordinates
    old_x = pb_atom[0]
    old_y = pb_atom[1]
    old_z = pb_atom[2]

    # distance in A. Margin determines which atoms will be reproduced by pbc
    margin = 4
    #determine low_boundary and high_boundary
    low_bound_x = margin
    high_bound_x = latpar[0] - margin
    low_bound_y = margin
    high_bound_y = latpar[1] - margin
    low_bound_z = margin
    high_bound_z = latpar[2] - margin
    #make list of coordinates
    coordinates_x = [old_x]
    coordinates_y = [old_y]
    coordinates_z = [old_z]
    if old_x < low_bound_x:
        new_x = old_x + latpar[0]
        coordinates_x.append(new_x)
    if old_y < low_bound_y:
        new_y = old_y + latpar[1]
        coordinates_y.append(new_y)
    if old_z < low_bound_z:
        new_z = old_z + latpar[2]
        coordinates_z.append(new_z)
    if old_x > high_bound_x:
        new_x = old_x - latpar[1]
        coordinates_x.append(new_x)
    if old_y > high_bound_y:
        new_y = old_y - latpar[1]
        coordinates_y.append(new_y)
    if old_z > high_bound_z:
        new_z = old_z - latpar[2]
        coordinates_z.append(new_z)

    if len(coordinates_x) == 1 and len(coordinates_y) == 1 and len(coordinates_z):
        return [pb_atom]
    else:
        for x in range(0, len(coordinates_x)):
            for y in range(0, len(coordinates_y)):
                for z in range(0, len(coordinates_z)):
                    twin_atom = [coordinates_x[x], coordinates_y[y], coordinates_z[z]]
                    replicas.append(twin_atom)
        return replicas

def get_cell(line):
    # split the file
    # Remove spaces from line
    line_without_spaces = line.strip()
    # Separate elements with a space
    split_lines = line_without_spaces.split()
    float_list = [float(x) for x in split_lines]
    return float_list


def get_coord(line):
    # split the file
    # Remove spaces from line
    line_without_spaces = line.strip()
    # Separate elements with a space
    split_lines = line_without_spaces.split()
    float_list = [float(x) for x in split_lines]
    coordinates = float_list[2:]
    return coordinates

def angle_vector(pb_i_triad):
    #subtract pb-i atoms
    vector_a = np.subtract(pb_i_triad[0], pb_i_triad[1])
    vector_b = np.subtract(pb_i_triad[0], pb_i_triad[2])
    # product between vectors
    inner = np.inner(vector_a, vector_b)
    norms = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
    #search principal axis between three atoms
    # it must be a subtraction between vectors
    angle_direction = np.subtract(vector_a, vector_b)
    if abs(angle_direction[0]) > abs(angle_direction[1]) and abs(angle_direction[0]) > abs(angle_direction[2]):
        direction = "x"
    elif abs(angle_direction[1]) > abs(angle_direction[0]) and abs(angle_direction[1]) > abs(angle_direction[2]):
        direction = "y"
    elif abs(angle_direction[2]) > abs(angle_direction[0]) and abs(angle_direction[2]) > abs(angle_direction[1]):
        direction = "z"
    else:
        direction = "error"
    # calculate angle
    cos = inner / norms
    rad = np.arccos(np.clip(cos, -1.0, 1.0))
    deg = np.rad2deg(rad)
    return deg, direction

def save_to_file(data_to_save, where_to_save):
    # Open the file
    with open(where_to_save, 'w') as output:
        # Iterate over each line
        for line in data_to_save:
            num_to_string = str(line)
            # Save line to file
            output.write(num_to_string + '\n')

def save_layers(data_to_save, where_to_save):
    # Open the file
    with open(where_to_save, 'w') as output:
        # iterate over each list
        for list in range(0, len(data_to_save)):
            # Iterate over each line
            output.write("# layer" + str(list + 1) + '\n')
            for line in data_to_save[list]:
                num_to_string = str(line)
                # Save line to file
                output.write(num_to_string + '\n')
        output.write('\n')


def histogram(angles):
    #bins_histogram = np.histogram(angles, bins=15, range=(150, 180))
    bins_seq = [142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162,
                163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180]
    data_bin = plt.hist(angles, bins=bins_seq, density=True, facecolor='r')
    # plt.bar(, align='center')
    x_bins = data_bin[1].tolist()
    y_bins = data_bin[0].tolist()
    histogram_dat_file = []
    for bin in range(0, len(y_bins)):
        bin_hist = int(x_bins[bin] + x_bins[bin + 1]) * 0.5, y_bins[bin]
        bin_line = str(bin_hist[0]) + "   " + str(bin_hist[1])
        histogram_dat_file.append(bin_line)
    #plt.xlabel('Angle (deg)')
    #plt.ylabel('Frequency')
    #plt.title('Angles')
    #plt.xlim(140, 180)
    #plt.ylim(0, 0.2)
    #plt.grid(True)
    #plt.savefig('hist_' + name + '.png')
    return histogram_dat_file

def angles(pb_i_list):
    #
    # First, the routine gets the Pb-I bonded atoms indexes, and calculates angles. In the next steps,
    # knowing the atoms' indexes, only update the angle values.
    #
    num_timesteps = len(pb_i_list)
    # Remember index of neighbouring pb-i atoms
    # master_list stores I-index, Pbs-index, bond direction axis, z-layer
    master_list = []
    # angles_list stores angles without distinction of direction and z-layer
    angles_list_total = []
    # angle_list_n stores angles in a certain direction, without distinction of the z-layer
    angle_list_x = []
    angle_list_y = []
    angle_list_z = []
    # angle_n_layer stores angles in a certain direction, making distinction of the z-layer
    angle_x_layer = []
    angle_y_layer = []
    angle_z_layer = []
    layer_list_x = []
    layer_list_y = []
    layer_list_z = []

    # i use this radius to separate layers and as a cutoff to select pb atoms close to I atoms
    radius_layer = 4
    latpar = []
    # For the first timestep, get first neighbours
    first_timestep_coordinates = pb_i_list[0]
    i_list = first_timestep_coordinates[1]
    pb_list = first_timestep_coordinates[0]
    for i in range(0, len(i_list)):
        bond_pbi_bond = []
        i_atom = i_list[i]
        i_z_value = i_atom[2]
        # select z layer
        z_layer = int(round(i_z_value, 0))
        i_index = int(i)
        atoms_indexes = [i_index]
        bond_pbi_bond.append(i_atom)
        for num_order_pb in range(0, len(pb_list)):
            #iterates over copies of Pb atom, stored by pbc
            pb_copies = pb_list[num_order_pb]
            for pb in pb_copies:
                pb_atom = pb
                difference_vector = np.subtract(i_atom, pb_atom)
                module_pb_i = np.linalg.norm(difference_vector)
                pb_index = int(num_order_pb)
                # now, I get the pb atoms next to the I one
                if module_pb_i < radius_layer:
                    #store pb-i indexes:
                    atoms_indexes.append(pb_index)
                    bond_pbi_bond.append(pb_atom)
                    # if i already have the I atom, and both Pb neighbours
                    if len(bond_pbi_bond) == 3:
                        # calculate angles
                        angle_between_atoms = angle_vector(bond_pbi_bond)
                        # separate atoms by axis
                        master_list_element = []
                        if angle_between_atoms[1] == "x":
                            # SAVE DIRECTION
                            # save angle in list of x angles
                            angle_list_x.append(angle_between_atoms[0])
                            # arrange in z-layers (I'm only interested in z layers)
                            # if the list is empty
                            if len(layer_list_x) == 0:
                                layer_list_x.append(z_layer)
                                angle_x_layer = [[z_layer, angle_between_atoms[0]]]
                                master_list_element = [atoms_indexes, z_layer, "x"]
                            # if the list is not empty, store the angle in a preexistent layer, or create a new one
                            else:
                                # counter to select preexistent layer or new
                                count = 0
                                for element in range(0, len(layer_list_x)):
                                    # check if the iodine atom is in preexistent layer
                                    layer_coord_z = layer_list_x[element]
                                    difx = abs(np.subtract(layer_coord_z, i_z_value))
                                    # if the layer exists, store the iodine atom
                                    if difx < radius_layer:
                                        master_list_element = [atoms_indexes, layer_coord_z, "x"]
                                        # workaround to avoid the creation of a new layer
                                        count = count + 1
                                        angle_x_layer[element].append(angle_between_atoms[0])
                                        break
                                    # if the layer doesn't exists, create a new one and store it
                                if count == 0:
                                    layer_list_x.append(z_layer)
                                    master_list_element = [atoms_indexes, z_layer, "x"]
                                    new_layer = [z_layer, angle_between_atoms[0]]
                                    angle_x_layer.append(new_layer)
                            master_list.append(master_list_element)
                        if angle_between_atoms[1] == "y":
                            # SAVE DIRECTION
                            # save angle in list of y angles
                            angle_list_y.append(angle_between_atoms[0])
                            # arrange in z-layers (I'm only interested in z layers)
                            # if the list is empty
                            if len(layer_list_y) == 0:
                                layer_list_y.append(z_layer)
                                angle_y_layer = [[z_layer, angle_between_atoms[0]]]
                                master_list_element = [atoms_indexes, z_layer, "y"]
                            # if the list is not empty, store the angle in a preexistent layer, or create a new one
                            else:
                                # counter to select preexistent layer or new
                                count = 0
                                for element in range(0, len(layer_list_y)):
                                    # check if the iodine atom is in preexistent layer
                                    layer_coord_z = layer_list_y[element]
                                    dify = abs(np.subtract(layer_coord_z, i_z_value))
                                    # if the layer exists, store the iodine atom
                                    if dify < radius_layer:
                                        master_list_element = [atoms_indexes, layer_coord_z, "y"]
                                        # workaround to avoid the creation of a new layer
                                        count = count + 1
                                        angle_y_layer[element].append(angle_between_atoms[0])
                                        break
                                    # if the layer doesn't exists, create a new one and store it
                                if count == 0:
                                    layer_list_y.append(z_layer)
                                    master_list_element = [atoms_indexes, z_layer, "y"]
                                    new_layer = [z_layer, angle_between_atoms[0]]
                                    angle_y_layer.append(new_layer)
                            master_list.append(master_list_element)
                        if angle_between_atoms[1] == "z":
                            # SAVE DIRECTION
                            # save angle in list of z angles
                            angle_list_z.append(angle_between_atoms[0])
                            # arrange in z-layers (I'm only interested in z layers)
                            # if the list is empty
                            if len(layer_list_z) == 0:
                                layer_list_z.append(z_layer)
                                angle_z_layer = [[z_layer, angle_between_atoms[0]]]
                                master_list_element = [atoms_indexes, z_layer, "z"]
                            # if the list is not empty, store the angle in a preexistent layer, or create a new one
                            else:
                                # counter to select preexistent layer or new
                                count = 0
                                for element in range(0, len(layer_list_z)):
                                    # check if the iodine atom is in preexistent layer
                                    layer_coord_z = layer_list_z[element]
                                    difz = abs(np.subtract(layer_coord_z, i_z_value))
                                    # if the layer exists, store the iodine atom
                                    if difz < radius_layer:
                                        master_list_element = [atoms_indexes, layer_coord_z, "z"]
                                        # workaround to avoid the creation of a new layer
                                        count = count + 1
                                        angle_z_layer[element].append(angle_between_atoms[0])
                                        break
                                    # if the layer doesn't exists, create a new one and store it
                                if count == 0:
                                    layer_list_z.append(z_layer)
                                    master_list_element = [atoms_indexes, z_layer, "z"]
                                    new_layer = [z_layer, angle_between_atoms[0]]
                                    angle_z_layer.append(new_layer)
                            master_list.append(master_list_element)
                        angles_list_total.append(angle_between_atoms[0])
    # now get angles from the remaining timesteps, using the memory_list_neighbour_pbi
    # I exclude the first timestep to avoid counting twice
    for timestep in pb_i_list[1:]:
        coordinates_timestep = timestep
        for trio in master_list:
            #from the master list, I get the atoms and calculate angles in each timestep
            i_atom = coordinates_timestep[1][trio[0][0]]
            z_value_i = round(i_atom[2], 0)
            # As pb atoms have instances due to pbc, I have to iterate on them
            pb_atom_1 = coordinates_timestep[0][trio[0][1]]
            for replicated_atom_1 in pb_atom_1:
                difference_vector = np.subtract(i_atom, replicated_atom_1)
                module_pb_i = np.linalg.norm(difference_vector)
                if module_pb_i < radius_layer:
                    pb_1 = replicated_atom_1
            pb_atom_2 = coordinates_timestep[0][trio[0][2]]
            for replicated_atom_2 in pb_atom_2:
                difference_vector = np.subtract(i_atom, replicated_atom_2)
                module_pb_i = np.linalg.norm(difference_vector)
                if module_pb_i < radius_layer:
                    pb_2 = replicated_atom_2
            three_atoms = [i_atom, pb_1, pb_2]
            angle_between_atoms = angle_vector(three_atoms)
            # store atoms by directions and planes
            angles_list_total.append(angle_between_atoms[0])
            if trio[2] == "x":
                #compares the layer of the atom with those stored from the first timestep
                angle_list_x.append(angle_between_atoms[0])
                for layer in range(0, len(angle_x_layer)):
                    layer_stored_x = angle_x_layer[layer][0]
                    # when the layers match, store the angle in the corresponding layer
                    difx = abs(layer_stored_x - z_value_i)
                    if difx < radius_layer:
                        angle_x_layer[layer].append(angle_between_atoms[0])
                        break
            if trio[2] == "y":
                #compares the layer of the atom with those stored from the first timestep
                angle_list_y.append(angle_between_atoms[0])
                for layer in range(0, len(angle_y_layer)):
                    layer_stored_y = angle_x_layer[layer][0]
                    # when the layers match, store the angle in the corresponding layer
                    dify = abs(layer_stored_y - z_value_i)
                    if dify < radius_layer:
                        angle_y_layer[layer].append(angle_between_atoms[0])
                        break
            if trio[2] == "z":
                #compares the layer of the atom with those stored from the first timestep
                angle_list_z.append(angle_between_atoms[0])
                for layer in range(0, len(angle_z_layer)):
                    layer_stored_z = angle_z_layer[layer][0]
                    # when the layers match, store the angle in the corresponding layer
                    difz = abs(layer_stored_z - z_value_i)
                    if difz < radius_layer:
                        angle_z_layer[layer].append(angle_between_atoms[0])
                        break
    # creates a list of angles per layer, counting from the top surface and including it in the first layer
    angles_per_layer = []
    for layer in range(0, len(angle_z_layer)):
        sum_layer = []
        sum_layer.extend(angle_x_layer[layer])
        sum_layer.extend(angle_y_layer[layer])
        sum_layer.extend(angle_z_layer[layer])
        angles_per_layer.append(sum_layer)
    # now, I store the opposite surface, which only have x-y angles
    sum_layer = []
    sum_layer.extend(angle_x_layer[-1])
    sum_layer.extend(angle_y_layer[-1])
    angles_per_layer.append(sum_layer)
    return angles_list_total, angle_list_x, angle_list_y, angle_list_z, angle_x_layer, angle_y_layer, angle_z_layer, angles_per_layer

def main(dlpoly_file):
    # read HISTORY file
    pb_iodide_atoms_time = read_file(dlpoly_file)
    # get angles form HISTORY file
    pb_i_angles = angles(pb_iodide_atoms_time)
    # make histogram from angles_total
    histogram_total = histogram(pb_i_angles[0])
    angle_file = "angle-total.dat"
    save_to_file(histogram_total, angle_file)
    histogram_x = histogram(pb_i_angles[1])
    angle_file_x = "angle-x.dat"
    save_to_file(histogram_x, angle_file_x)
    histogram_y = histogram(pb_i_angles[2])
    angle_file_y = "angle-y.dat"
    save_to_file(histogram_y, angle_file_y)
    histogram_z = histogram(pb_i_angles[3])
    angle_file_z = "angle-z.dat"
    save_to_file(histogram_z, angle_file_z)
    # layers_x
    angle_histogram_x = []
    for m in range(0, len(pb_i_angles[4])):
        layer = pb_i_angles[4][m]
        histogram_layer = histogram(layer)
        angle_histogram_x.append(histogram_layer)
    angle_layers_x_file = "angles_layers_x.dat"
    save_layers(angle_histogram_x, angle_layers_x_file)
    # layers y
    angle_histogram_y = []
    for n in range(0, len(pb_i_angles[5])):
        layer = pb_i_angles[5][n]
        histogram_layer = histogram(layer)
        angle_histogram_y.append(histogram_layer)
    angle_layers_y_file = "angles_layers_y.dat"
    save_layers(angle_histogram_y, angle_layers_y_file)
    # layers y
    angle_histogram_z = []
    for o in range(0, len(pb_i_angles[6])):
        layer = pb_i_angles[6][o]
        histogram_layer = histogram(layer)
        angle_histogram_z.append(histogram_layer)
    angle_layers_z_file = "angles_layers_z.dat"
    save_layers(angle_histogram_z, angle_layers_z_file)
    # save angles_per_layer
    angle_layer_hist = []
    for p in range(0, len(pb_i_angles[7])):
        layer = pb_i_angles[7][p]
        histogram_layer = histogram(layer)
        angle_layer_hist.append(histogram_layer)
    angle_layers_file = "angles_by_layers.dat"
    save_layers(angle_layer_hist, angle_layers_file)




dlpoly_file = './phas2.dat'
# Always import python libraries "at the top" of the program (physically, at the "end" of it)
# bc it is available to any subroutine.
# numpy bc I need to make some calculations
import numpy as np
# matplotlib to plot histograms
import matplotlib.pyplot as plt
z_angle = main(dlpoly_file)
