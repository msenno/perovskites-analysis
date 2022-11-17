#modified from ma_slices_history, now takes y-slices
#Terminar de arreglar
def read_file(file_path):
    # read the file
    with open(file_path) as f:
        # Get all file lines
        lines = f.readlines()
    return lines


def get_coord(line):
    # split the file
    # Remove spaces from line
    line_without_spaces = line.strip()
    # Separate elements with a space
    split_lines = line_without_spaces.split()
    float_list = [float(x) for x in split_lines]
    return float_list

def projection(c_coord, n_coord, latpar):
    # corrects the pbc in molecules
    # This tiny step takes care of the Periodic Boundary Conditions (pbc). When a C or N atom
    # appears in the opposite side of the cell, this subroutine corrects its position.
    # limit_boundary select the minimum value in each coord. Below that limit, atoms will be rearranged
    # to the opposite side of the cell. It moves the atoms from one face of the cell to another.

    v_x = abs(float(n_coord[0]) - float(c_coord[0]))
    v_y = abs(float(n_coord[1]) - float(c_coord[1]))
    v_z = abs(float(n_coord[2]) - float(c_coord[2]))

    if v_x > 0.5*latpar[0]:
        if c_coord[0] < n_coord[0]:
            c_coord[0] = c_coord[0] + latpar[0]
        else:
            n_coord[0] = n_coord[0] + latpar[0]
        v_x = abs(float(n_coord[0]) - float(c_coord[0]))

    if v_y > 0.5*latpar[1]:
        if c_coord[1] < n_coord[1]:
            c_coord[1] = c_coord[1] + latpar[1]
        else:
            n_coord[1] = n_coord[1] + latpar[1]
        v_y = abs(float(n_coord[1]) - float(c_coord[1]))
    # get module of vector
    if 0.5*c_coord[0] + 0.5*n_coord[0] < 4:
        c_coord[0] = c_coord[0] + latpar[0]
        n_coord[0] = n_coord[0] + latpar[0]
    if 0.5*c_coord[1] + 0.5*n_coord[1] < 4:
        c_coord[1] = c_coord[1] + latpar[1]
        n_coord[1] = n_coord[1] + latpar[1]

    position_molecule_x = (c_coord[0] + n_coord[0])/2
    position_molecule_y = (c_coord[1] + n_coord[1])/2
    position_molecule_z = (c_coord[2] + n_coord[2])/2

    module3_d = (v_x ** 2 + v_y ** 2 + v_z ** 2)**0.5
    cosine = round((v_z / module3_d), 4)
    molecule_data = [position_molecule_x, position_molecule_z, v_x, v_z, cosine, position_molecule_y]
    return molecule_data


def extract_data(line):
    # Lee las coordenadas de C y N
    lines = read_file(line)
    # Parameters of file
    num_lines = len(lines)  # Number of lines in a timestep
    # Create list of coordinates
    all_molecules_data = []
    timestep_list = []
    layer_list = []
    radius = 2.5
    latpar = []
    # For loop, takes every C and N and saves their coordinates
    for i in range(0, num_lines):
        # I need to save data before the program reaches de last line, so I save here
        if i == num_lines - 1:
            timestep_list.append(layer_list)
        # Goes line by line
        textline = lines[i].strip().split()
        # Use this line with HISTORY files
        # before restart temporal_list, save data in timestep_list
        if "TIMESTEP" in textline:
            #read lattice parameters in case of rearrange by pbc
            x_lattice = get_coord(lines[i + 5])[1]
            y_lattice = get_coord(lines[i + 6])[1]
            z_lattice = get_coord(lines[i + 7])[1] - get_coord(lines[i + 7])[0]
            latpar = [float(x_lattice), float(y_lattice), float(z_lattice)]
            if len(layer_list) > 0:
                timestep_list.append(layer_list)
            # reboot temporal list
            layer_list = []
        if len(textline) > 3:
            if textline[1] == "3":
                c_line = lines[i]
                c_coord = get_coord(c_line)[2:]
                n_line = lines[i + 1]
                n_coord = get_coord(n_line)[2:]
                #projections take the coordinates of n and c and returns positions, vectors and cos(theta)
                #pos_x pos_z v_x v_z cos_theta pos_y
                projections = projection(c_coord, n_coord, latpar)
                # select x layer
                y_coordinate = round(projections[5], 1)
                if len(layer_list) == 0:  # turnaround when layer_list is empty
                    layer_list.append([y_coordinate])  # add first value
                    layer_list[0].append(projections)
                else:  # add or select layer
                    count = 0
                    for element in range(0, len(layer_list)):
                        # if the absolute difference between z coordinate of molecule and every
                        # element of the list is lower than a cutoff value, I place it in one
                        # layer, if not, I place it in a new layer.
                        difference = abs(y_coordinate - layer_list[element][0])  # make difference
                        if difference < radius:  # if lower than radius, add to same layer
                            layer_list[element].append(projections)
                            count = count + 1
                    if count == 0:
                        layer_list.append([y_coordinate])
                        # It is a new element, and its index value is len(list) - 1. If the list has 3 elements,
                        # the new element has been indexed as 2 (remember, python reads 0,1,2)
                        layer_list[len(layer_list) - 1].append(projections)
                    else:
                        continue
    return timestep_list


def average(timesteps):
    # define average in time
    num_dipoles = len(timesteps[0][0])
    averages_list_final = []
    len_step = len(timesteps[0])
    # Factor: enlarges orientation vector
    factor = 3
    for layer in range(0, len_step):
        layer_list = []
        for dipole in range(1, num_dipoles):
            z_projection = 0
            pos_molx = 0
            pos_molz = 0
            vec_x = 0
            vec_z = 0
            for tstep in range(0, len(timesteps)):
                projection_z = timesteps[tstep][layer][dipole][4] / len(timesteps)
                coord_pos_mol_x = timesteps[tstep][layer][dipole][0] / len(timesteps)
                coord_pos_mol_z = timesteps[tstep][layer][dipole][1] / len(timesteps)
                vec_cation_x = timesteps[tstep][layer][dipole][2] / len(timesteps) * factor
                vec_cation_z = timesteps[tstep][layer][dipole][3] / len(timesteps) * factor
                z_projection = z_projection + projection_z
                pos_molx = pos_molx + coord_pos_mol_x
                pos_molz = pos_molz + coord_pos_mol_z
                vec_x = vec_x + vec_cation_x
                vec_z = vec_z + vec_cation_z
            position_vectors_proj = [pos_molx, pos_molz, vec_x, vec_z, z_projection]
            layer_list.append(position_vectors_proj)
        averages_list_final.append(layer_list)
    return averages_list_final

def save_to_file(data_to_save, where_to_save):
    # Open the file
    # Iterate over each line
    for line in range(0, len(data_to_save)):
        # Save line to file
        subfile = where_to_save + '-layer' + str('{:0>3}'.format(line + 1)) + '.dat'
        with open(subfile, 'w') as output:
            for element in range(1, len(data_to_save[line])):
                string_to_print = '\t'.join(str(elem) for elem in data_to_save[line][element])
                output.write(string_to_print + '\n')


def save_average(data_to_save, where_to_save):
    # Open the file
    # Iterate over each line
    with open(where_to_save, 'w') as output:
        for element in range(0, len(data_to_save)):
             # Save line to file
            line_to_save = data_to_save[element]
            string_line = '\t'.join(str(round(elem, 5)) for elem in line_to_save)
            output.write(string_line + '\n')
        output.write('\n')

def main(file_path):
    print("Choose which type of time processing you want")
    #time_process = input("1) Average     2) By timestep \n")
    time_process = str(1)
    timesteps = extract_data(file_path)
    if time_process == str(1):
        aver_time = average(timesteps)
        for j in range(0, len(aver_time)):
            file_to_save_result = './aver-ylayer' + str('{:0>2}'.format(j + 1)) + '.dat'
            save_average(aver_time[j], file_to_save_result)
    elif time_process == str(2):
        # Project vector in an axis
        # Save data to file
        for j in range(0, len(timesteps)):
            file_to_save_result = './tstep' + str('{:0>3}'.format(j + 1))
            save_to_file(timesteps[j], file_to_save_result)


layer_folder_path = ['phasetransition2.dat']
for folder_path in layer_folder_path:
    # open program "cos_theta_phi", which calculates angles and cosine
    main(folder_path)
