#!/usr/bin/env python
import numpy as np
import time 
import pdb_read.math_p as mp

def get_normal_distribution(sampleNo):
    sigma = 5
    np.random.seed(int(time.time()))
    s = np.random.normal(0, sigma, sampleNo)
    return s
def get_sidechain_dhddict():
    filename = "file/res.txt"
    with open (filename) as f:
        lines = f.read().strip()
        list_res = lines.split("\n\n")
        res_nums = len(list_res)
        dhd_atom_dict = {}
        for i in range(res_nums):
            dhd_info = list_res[i].split("\n")
            dhd_num = int(dhd_info[3][4])
            res_name = dhd_info[0][-3:]
            dhd_atoms = dhd_info[-dhd_num:]
            for j,atom in enumerate(dhd_atoms):
                atom = str(atom).split()
                dhd_atoms[j] = atom        
            #print(dhd_atoms)
            dhd_atom_dict[res_name] = dhd_atoms
    return dhd_atom_dict

def new_side_res_norm_dhd(rot_res, dhd_dict, res_rotnum_dict, res_atomnum_dict):
    #dict_res_dhd = {}
    aa_num = len(rot_res)
    for n in range(aa_num):
        res_name1 = rot_res[n][0][0][3]
        #res对应的侧链dhd列表
        rot_atom_list = dhd_dict[res_name1]
        #res的rot数
        rot_num = res_rotnum_dict[res_name1]
        #res包含的atom数
        res_atom_num = res_atomnum_dict[res_name1]
        #print(res_name1)
        if res_name1 == 'GLY' or res_name1 == 'ALA':
            continue
        elif res_name1 == 'ARG':
            for i in range(rot_num): 
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg_cood = rot_res[n][i][3][-1]
                cd_cood = rot_res[n][i][4][-1]
                ne_cood = rot_res[n][i][5][-1]
                cz_cood = rot_res[n][i][6][-1]                
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                cg_cd = cg_cood - cd_cood
                cd_ne = cd_cood - ne_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg, cg_cd, cd_ne]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2,3,4]))
                for x1 in range(3,7):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(4,7):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])
                for x3 in range(5,7):
                    cood = rot_res[n][i][x3][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[2], sheta_list[2])
                for x4 in range(6,7):
                    cood = rot_res[n][i][x4][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[3], sheta_list[3])
        elif res_name1 == 'ASN' or res_name1 == 'ASP':
            for i in range(rot_num):
                n1 = rot_res[n][i] 
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg_cood = rot_res[n][i][3][-1]
                od_cood = rot_res[n][i][4][-1]                
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2]))
                for x1 in range(3,5):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(4,5):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])
        elif res_name1 == 'CYS':
            for i in range(rot_num):
                n1 = rot_res[n][i] 
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                sg_cood = rot_res[n][i][3][-1]                
                ca_cb = ca_cood - cb_cood
                fa = ca_cb/np.linalg.norm(ca_cb)
                sheta = get_normal_distribution(10)[1]
                cood = rot_res[n][i][3][-1]
                rot_res[n][i][x1][-1] = mp.rotation(cood, fa, sheta)

        elif res_name1 == 'GLN' or res_name1 == 'GLU':
            for i in range(rot_num):
                n1 = rot_res[n][i] 
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg_cood = rot_res[n][i][3][-1]
                cd_cood = rot_res[n][i][4][-1]
                oe_cood = rot_res[n][i][5][-1]
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                cg_cd = cg_cood - cd_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg, cg_cd]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2,3]))
                for x1 in range(3,6):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(4,6):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])
                for x3 in range(5,6):
                    cood = rot_res[n][i][x3][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[2], sheta_list[2])
        elif res_name1 == 'HIS' or res_name1 == 'HSD' or res_name1 == 'HSE' or res_name1 == 'HSP':
            for i in range(rot_num):
                n1 = rot_res[n][i] 
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cd_cood = rot_res[n][i][3][-1]
                cg_cood = rot_res[n][i][4][-1]
                ne_cood = rot_res[n][i][5][-1]
                nd_cood = rot_res[n][i][6][-1]                
                ce_cood = rot_res[n][i][7][-1]              
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2]))
                for x1 in range(3,8):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(4,8):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])
        elif res_name1 == 'ILE':
            for i in range(rot_num):
                n1 = rot_res[n][i] 
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg_cood = rot_res[n][i][3][-1]
                cd_cood = rot_res[n][i][4][-1]                
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2]))
                for x1 in range(3,5):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(4,5):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])

        elif res_name1 == 'LYS':
            for i in range(rot_num):
                n1 = rot_res[n][i] 
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg_cood = rot_res[n][i][3][-1]
                cd_cood = rot_res[n][i][4][-1]
                ce_cood = rot_res[n][i][5][-1]
                nz_cood = rot_res[n][i][6][-1]                
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                cg_ce = cg_cood - ce_cood
                cd_nz = cd_cood - nz_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg, cg_cd, cd_ne]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2,3,4]))
                for x1 in range(3,7):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(4,7):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])
                for x3 in range(5,7):
                    cood = rot_res[n][i][x3][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[2], sheta_list[2])
                for x4 in range(6,7):
                    cood = rot_res[n][i][x4][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[3], sheta_list[3])
        elif res_name1 == 'MET':
            for i in range(rot_num):
                n1 = rot_res[n][i] 
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg_cood = rot_res[n][i][3][-1]
                sd_cood = rot_res[n][i][4][-1]
                ce_cood = rot_res[n][i][5][-1]               
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                cg_sd = cg_cood - cd_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg, cg_sd]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2,3]))
                for x1 in range(3,6):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(4,6):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])
                for x3 in range(5,6):
                    cood = rot_res[n][i][x3][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[2], sheta_list[2])
        elif res_name1 == 'PHE':
            for i in range(rot_num):
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg_cood = rot_res[n][i][3][-1]
                cd_cood = rot_res[n][i][4][-1]               
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg, cg_cd, cd_ne]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2,3,4]))
                for x1 in range(3,9):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(4,9):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x2][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])
        elif res_name1 == 'PRO':
            for i in range(rot_num):
                n_cood = rot_res[n][i][0][-1]
                cd_cood = rot_res[n][i][1][-1]
                ca_cood = rot_res[n][i][2][-1]
                cb_cood = rot_res[n][i][3][-1]
                cg_cood = rot_res[n][i][4][-1]               
                ca_cb = ca_cood - cb_cood
                fa = ca_cb/np.linalg.norm(ca_cb)
                sheta = get_normal_distribution(10)[1]
                cood = rot_res[n][i][4][-1]
                rot_res[n][i][4][-1] = mp.rotation(cood, fa, sheta)
        elif res_name1 == 'SER':
            for i in range(rot_num):
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                og_cood = rot_res[n][i][3][-1]               
                ca_cb = ca_cood - cb_cood
                fa = ca_cb/np.linalg.norm(ca_cb)
                sheta = get_normal_distribution(10)[1]
                cood = rot_res[n][i][3][-1]
                rot_res[n][i][3][-1] = mp.rotation(cood, fa, sheta)
        elif res_name1 == 'THR':
            for i in range(rot_num):
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                og_cood = rot_res[n][i][3][-1]               
                ca_cb = ca_cood - cb_cood
                fa = ca_cb/np.linalg.norm(ca_cb)
                sheta = get_normal_distribution(10)[1]
                for x1 in range(3,5):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood, fa, sheta)
        elif res_name1 == 'TRP':
            for i in range(rot_num):
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg_cood = rot_res[n][i][3][-1]
                cd2_cood = rot_res[n][i][7][-1]               
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2]))
                for x1 in range(3,11):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(7,11):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x2][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])
        elif res_name1 == 'TYR':
            for i in range(rot_num):
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg_cood = rot_res[n][i][3][-1]
                cd1_cood = rot_res[n][i][4][-1]               
                ca_cb = ca_cood - cb_cood
                cb_cg = cb_cood - cg_cood
                fa_list = list(map(lambda x: x/np.linalg.norm(x), [ca_cb, cb_cg]))
                sheta_list = list(map(lambda x: get_normal_distribution(10)[x], [1,2]))
                for x1 in range(3,9):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood,fa_list[0], sheta_list[0])
                for x2 in range(4,9):
                    cood = rot_res[n][i][x2][-1]
                    rot_res[n][i][x2][-1] = mp.rotation(cood,fa_list[1], sheta_list[1])
        elif res_name1 == 'VAL':
            for i in range(rot_num):
                n_cood = rot_res[n][i][0][-1]
                ca_cood = rot_res[n][i][1][-1]
                cb_cood = rot_res[n][i][2][-1]
                cg1_cood = rot_res[n][i][3][-1]               
                ca_cb = ca_cood - cb_cood
                fa = ca_cb/np.linalg.norm(ca_cb)
                sheta = get_normal_distribution(10)[1]
                for x1 in range(3,5):
                    cood = rot_res[n][i][x1][-1]
                    rot_res[n][i][x1][-1] = mp.rotation(cood, fa, sheta)
    return rot_res





















