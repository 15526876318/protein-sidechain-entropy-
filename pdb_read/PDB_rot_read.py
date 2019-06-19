#!/usr/bin/env python
import numpy as np
res_names = (
    'ALA',
    'ARG',
    'ASN',
    'ASP',
    'CYS',
    'GLN',
    'GLU',
    'GLY',
    'HIS',
    'ILE',
    'LEU',
    'LYS',
    'MET',
    'PHE',
    'PRO',
    'SER',
    'THR',
    'TRP',
    'TYR',
    'VAL')

class PDB_rot_read(object):
    ''' read rotamer database'''
    def __init__(self):
        self.rot_tuple = []
        self.resrotnum = np.zeros((20), dtype=np.int)
        self.resatmnum = np.zeros((20), dtype=np.int)
        self.res_atomnum_dict = {}
        self.res_rotnum_dict = {}
        self.main_res_dict = {}
        self.main_res_num = []
        self.main_res_name = []
        self.main_res_atomnum = []
        self.main_res_rotnum = []
        self.rot_res = [[] for num in range(20)]

    def read_pdbfile(self, file_path):
        append = self.rot_tuple.append
        with open(file_path, "r") as pdb_file:
            for line_str in pdb_file:
                line_str = line_str.strip()
                line_str = ' '.join(line_str.split())
                id_name, serial_num, atom_name, res_name, res_num, x, y, z, bfactor, num, st = line_str.split(' ')
                coord = np.array((x, y, z), dtype=np.float)
                # print(atom_name[0])
                if atom_name[0] != 'H':
                    pdb_tuple = [
                        int(serial_num),
                        atom_name,
                        int(res_num),
                        res_name,
                        coord]
                    append(pdb_tuple)

    def get_information_rot(self, res_names, norm = None):
        rot_serials = len(self.rot_tuple)
        resnum = 0
        n = 0
        line = self.rot_tuple
        for i in range(rot_serials - 1):
            if line[i][3] != line[i + 1][3]:
                self.resrotnum[resnum] = line[i][2]
                resnum += 1
        self.resrotnum[19] = line[rot_serials - 1][2]

        for num in range(20):
            atom = 0
            for i in range(rot_serials):
                if line[i][3] == res_names[num]:
                    atom += 1
            self.resatmnum[num] = atom / self.resrotnum[num]
        for i in range(20):
            self.rot_res[i] = [[] for num1 in range(self.resrotnum[i])]
        for x in range(20):
            for y in range(self.resrotnum[x]):
                for z in range(self.resatmnum[x]):
                    self.rot_res[x][y].append(line[n])
                    # print(self.rot_res[x][y][z])
                    n += 1
        for i in range(20):
            self.res_atomnum_dict[res_names[i]] = self.resatmnum[i]
            self.res_rotnum_dict[res_names[i]] = self.resrotnum[i]
        
        if norm:
            dhd_dict = get_sidechain_dhddict()
            new_side_res_norm_dhd(self.rot_res, dhd_dict, self.res_rotnum_dict, self.res_atomnum_dict)


    def mainchain_group(self, res_num, mainchain, res_names, start_resnum):
        #self.check_res(mainchain)
        main_res_group = [[] for num in range(res_num)]
        for j in range(res_num):
            for i in range(len(mainchain)):
                if int(mainchain[i][2]) == j + start_resnum + 1:
                    main_res_group[j].append(mainchain[i]) ###############
        self.main_res_group = main_res_group
        #try:
            #if len(main_res_group) !=
        #resname11 = main_res_group[i][0][3]  # list
        for i in range(res_num):
                self.main_res_name.append(main_res_group[i][0][3])###
        for m,res in enumerate(res_names):
            self.main_res_dict[res] = m 
        self.main_res_num = list(map(lambda m: self.main_res_dict[m], self.main_res_name))
        self.main_res_atomnum = list(map(lambda m: self.res_atomnum_dict[m], self.main_res_name))
        self.main_res_rotnum = list(map(lambda m: self.res_rotnum_dict[m], self.main_res_name))
    def check_res(self, mainchain):
        "check the residue whether lack cood"
        check = False
        for i,atom_info in enumerate(mainchain):
            if i > 0 and atom_info[2] - mainchain[i-1][2]> 1:
                print(atom_info[2], mainchain[i-1][2],i)
                check = True
        if check:
            exit()
