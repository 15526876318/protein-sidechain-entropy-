#!/usr/bin/env python
import numpy as np
from numpy import reshape, tile
import math
import time
import sys
import warnings
import copy
import gzip
import re
import os
import pdb_read.math_p as mp
from pdb_read.PDB_rot_read import PDB_rot_read
from pdb_read.math_p import L_MO_ab
from pdb_read.PDBExceptions import PDBConstructionException
from pdb_read import SCE93_caleach2
import pdb_read.SCE10 as sc
from side_norm_dis import get_sidechain_dhddict, new_side_res_norm_dhd


vdw = {"H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.85}
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

class Fit_it(object):
    def __init__(
            self,
            res_num,
            main_res_name,
            main_res_group,
            rot_res,
            resrotnum,
            resatmnum,
            res_names,
            res_rotnum_dict,
            res_atomnum_dict,
            vdw,
            w):
        self.res_num = res_num
        self.main_res_group = main_res_group
        self.rot_res = rot_res
        self.resrotnum = resrotnum
        self.resatmnum = resatmnum
        self.res_names = res_names
        self.res_rotnum_dict = res_rotnum_dict
        self.res_atomnum_dict = res_atomnum_dict
        self.vdw = vdw
        self.main_res_name = main_res_name
        #self.atomnumc = np.zeros((res_num, 35), dtype = np.int)

    def res_fit(self):
        res_rot_fited = [[] for num in range(self.res_num)]
        self.atomnumc = np.zeros((self.res_num,34),dtype= np.int32)
        rot_res = copy.deepcopy(self.rot_res)                
        for i in range(self.res_num):
            for j in range(20):
                if self.main_res_group[i][0][3] == self.res_names[j]:   ###############################
                    res_rot_fited[i] = [[] for num2 in range(self.resrotnum[j])]
                    for k in range(self.resrotnum[j]):
                        res_rot_fited[i][k] = rot_res[j][k]
        self.res_rot_fited = res_rot_fited
        for i in range(self.res_num):
            self.mainchain_coord = self.main_res_group[i]
            resname = self.main_res_name[i]
            rotnum = self.res_rotnum_dict[resname]
            # print(rotnum)
            for k in range(rotnum):
                self.sidechain_coord = copy.deepcopy(self.res_rot_fited[i][k])
                if resname == 'GLY':
                    self.res_rot_fited[i][k] = self.main_res_group[i]
                else:
                    self.res_rot_fited[i][k] = self._fit_angle()
                    #self.atomnumc[i][k] = res_atomnum_dict[self.res_rot_fited[i][k][0][3]]
                    self.atomnumc[i][k] = self.res_atomnum_dict[self.res_rot_fited[i][k][0][3]]
    def _fit_angle(self):
        bn = self.sidechain_coord #dedao rotamer de zuobiao(bn)
        ar = self.mainchain_coord #pdb wenjian zhong de zuobiao(ar)
        #得到主链坐标
        try:
            ar_n = ar[0][4]
            ar_ca = ar[1][4]
            ar_c = ar[2][4]
            ar_o = ar[3][4]
            ar_cb = ar[4][4]
        except IndexError:
            print("%s res %d loss some mainchain cood! " % (protein_name, ar[0][2]))
            exit()
        for j, bn1 in enumerate(bn):
            if bn1[1] == 'CA':
                bn_ca = bn1[4]
            elif bn1[1] == 'CB':
                bn_cb = bn1[4]
        bn_cb_ca = bn_cb - bn_ca
        ar_cb_ca = ar_cb - ar_ca
        #求得主链和侧链ca_cb键角
        sheta = mp.get_angle(ar_cb_ca, bn_cb_ca)
        #print(sheta)
        v1 = mp.L_cxyz(bn_cb_ca, ar_cb_ca)
        v1 = v1 / np.linalg.norm(v1)
        for j, bn1 in enumerate(bn):  #an faxiangliang weizou angle zhuandong
            bn[j][4] = mp.rotation(bn[j][4], v1, sheta)
            if bn[j][1] == 'N':
                bn_n = bn[j][4]
            elif bn[j][1] == 'CA':
                bn_ca = bn[j][4]
            elif bn[j][1] == 'CB':
                bn_cb = bn[j][4]
        #print(bn)
        b1 = bn_cb - bn_ca
        a2 = ar_ca - ar_n
        b2 = bn_n - bn_ca
        sheta1 = mp.get_dhd(a2, b1, b2)
        m_b1 = np.linalg.norm(b1) #b1 zhuanhua cheng danwei xiangliang 
        b1 = b1 / m_b1 
        for j, bn1 in enumerate(bn): #rao b1 zhuandong ermianjiao
            bn[j][4] = mp.rotation(bn[j][4], b1, -sheta1)
            bn[j][4] = bn[j][4] + ar_ca
            if bn[j][1] == 'N':
                bn[j][4] = ar_n
            elif bn[j][1] == 'CA':
                bn[j][4] = ar_ca
            elif bn[j][1] == 'C':
                bn[j][4] = ar_c
            elif bn[j][1] == 'CB':
                bn[j][4] = ar_cb
            elif bn[j][1] == 'O':
                bn[j][4] = ar_o            
        return bn

    def get_rc(self, res_atomnum_dict):
        res_rot_fited = self.res_rot_fited
        res_num = self.res_num
        res_rotnum_dict = self.res_rotnum_dict
        vdw = self.vdw
        mrc = np.zeros((res_num, 35), dtype=np.float)  # 侧链球心半径
        mrc_axyz = np.zeros((res_num, 35, 3), dtype=np.float)  # 　侧链球心坐标
        mrc_d = np.zeros((res_num, 35, 14, 3), dtype=np.float)  # 侧链原子坐标
        atom_rc = np.zeros((res_num, 14), dtype=np.float)  # 侧链的原子banjin
        for res in range(res_num):
            res_atomnum = res_atomnum_dict[res_rot_fited[res][0][0][3]]
            if res_rot_fited[res][0][0][3] == "GLY":
                continue
            res_atomnum_c = res_atomnum - 4
            for rot in range(res_rotnum_dict[self.res_rot_fited[res][0][0][3]]):
                if res_rot_fited[res][0][0][3] == "PRO":
                    atom_n_pro = (1, 3, 4)
                    xyz = [res_rot_fited[res][rot][atom][4] for atom in atom_n_pro]
                else:
                    xyz = [res_rot_fited[res][rot][atom][4] for atom in range(2, res_atomnum - 2)]
                xyz = np.array(xyz, dtype=np.float)
                for atom in range(res_atomnum_c):
                    mrc_d[res][rot][atom] = xyz[atom]
                distance = np.median(xyz, axis=0)  # celianqiuxin juli
                tr = [np.linalg.norm(distance - xyz[l]) for l in range(res_atomnum_c)]
                trm = max(tr)
                mrc_axyz[res][rot] = distance
                mrc[res][rot] = trm
            if res_rot_fited[res][0][0][3] == "PRO":
                n_pro = (1, 3, 4)
                atom = [res_rot_fited[res][0][atom][1][0] for atom in n_pro]
            else:
                atom = [res_rot_fited[res][0][atom][1][0] for atom in range(2, res_atomnum - 2)]
            atom1 = list(map(lambda m: vdw[m], atom))
            for j in range(len(atom1)):
                atom_rc[res][j] = atom1[j]
        self.mrc = mrc
        self.mrc_axyz = mrc_axyz
        self.mrc_d = mrc_d
        self.atom_rc = atom_rc

    def get_rm(self, res_atomnum_dict):
        res_rot_fited = self.res_rot_fited
        res_num = self.res_num
        vdw = self.vdw
        mr = np.zeros((res_num), dtype=np.float)
        mr_axyz = np.zeros((res_num, 3), dtype=np.float)
        mr_d = np.zeros((res_num, 4, 3), dtype=np.float)
        atom_r = np.zeros((res_num, 4), dtype=np.float)
        for res in range(res_num):
            n = res_atomnum_dict[res_rot_fited[res][0][0][3]]  # 不同氨基酸包含的原子数
            n1 = [0, 1, n - 2,-1]
            n2 = [0, 2, n - 2,-1]
            if self.main_res_name[res] == "PRO":
                xyz = [res_rot_fited[res][0][atom][4] for atom in n2]  # 主链的模拟半径
            else:
                xyz = [res_rot_fited[res][0][atom][4] for atom in n1]
            xyz = np.array(xyz, dtype=np.float)
            a = np.median(xyz, axis=0)
            tr = [np.linalg.norm(a - xyz[l]) for l in range(4)]
            trm = max(tr)
            mr_d[res] = xyz
            mr[res] = trm
            mr_axyz[res] = a
            if self.main_res_name[res] == "PRO":
                atom = [res_rot_fited[res][0][atom][1][0] for atom in n2]
            else:
                atom = [res_rot_fited[res][0][atom][1][0] for atom in n1]
            atom1 = list(map(lambda m: vdw[m], atom))
            #print(atom1)
            for j in range(len(atom1)):
                atom_r[res][j] = atom1[j]
        self.mr = mr
        self.mr_axyz = mr_axyz
        self.mr_d = mr_d
        self.atom_r = atom_r

def mc_rot(
        w,
        cp,
        res_num,
        atomnumc,
        main_res_atomnum,
        main_res_rotnum,
        mrc_axyz,
        mrc_d,
        mr_d,
        atom_r,
        atom_rc,
        mr,
        mrc,
        mainatom_r,
        mainatom_xyz,
        main_res_name,
        main_res_num,
        success_sum_array,
        success_num_array):
    np.random.seed(int(time.time()))
    cy_s_cs = np.zeros((cp, res_num), dtype=np.int32)
    wt = np.zeros((cp), dtype=np.float)
    rs = np.zeros((cp), dtype=np.int32)
    rm = np.zeros((cp), dtype=np.int32)
    #dhd_dict = get_sidechain_dhddict()##################################################
    #threads = []
    for nn in range(res_num):
        rsn = 0
        rmn = 0
        atom_num = main_res_atomnum[nn]
        rot_num = main_res_rotnum[nn]
        for j in range(cp):            
            #new_side_res_norm_dhd(self.rot_res, dhd_dict, res_rotnum_dict, self.res_atomnum_dict) 
            success_sum = 0
            for rot in range(rot_num):
                a = mrc_axyz[nn][rot]
                satom_xyz = mrc_d[nn][rot]  # celian de yaunzi zuobiao jihe
                mrcnk = mrc[nn][rot] 
                t = SCE93_caleach2.caleach(w, res_num, nn, atom_num, rot, j, atomnumc, mrc_axyz, 
                                            mrc_d, mrc, mr, mr_d, mrcnk, mainatom_xyz, atom_r, atom_rc, a, cy_s_cs, main_res_num)
                if t:
                    success_num_array[nn][success_sum] = rot
                    success_sum += 1
            success_sum_array[nn] = success_sum
            t1 = success_sum_array[nn]  # ke shengcheng de celiangouxiangshu
            if not t1:
                wt[j] = 0.0
                rm[rmn] = j
                rmn += 1  # sianjisuan geshu
            else:
                wt[j] = wt[j] + math.log(t1)
                rdi = np.random.randint(1, t1 + 1)  # gouxiangshu suoyin
                t2 = success_num_array[nn][rdi - 1]  # mc moni de dijizhong gouxiang
                cy_s_cs[j][nn] = t2
                rs[rsn] = j
                rsn += 1
        if not rsn:
            print("cannot generate %s" % nn)
            #return 0
        for j in range(rmn):
            try:
                t1 = np.random.randint(1, rsn + 1)
            except ValueError:
                #print("cannot genetate! rsn =%d" % rsn )
                exit()
            t2 = rs[t1]
            t3 = rm[j]
            cy_s_cs[t3] = cy_s_cs[t2]
            wt[t3] = wt[t2]
    wtt = 0.0
    for i in range(cp):
        wtt += wt[i]
    result = wtt / cp
    return result
    #print("%.2f" % result)

def Biopdb_read(filename):
    p = sc.PDBP_read(PERMISSIVE=True)
    s = p.get_structure(filename[-11:-7], filename)
    mainchain_all = []
    for c in s[0]:
        chain_id = c.id
        pdb_array =[]
        for i in c.get_atoms():
            res = i.get_parent()
            res_name= res.resname            
            res_num = res.id[1]
            info = [i.serial_number, i.fullname.strip(), res_num, res_name, i.coord, chain_id]
            pdb_array.append(info)
        mainchain1 = [
            row for row in pdb_array if row[1] in (
                'CA', 'N', 'C', 'CB', 'O')]
        #print(len(mainchain1),"1")
        p = len(mainchain1) - 1
        mainchain = [mainchain1[j] for j in range(
            p) if mainchain1[j][1] != mainchain1[j - 1][1]]
        mainchain.append(mainchain1[-1])
        #print(len(mainchain))
        mainchain_all.append(mainchain)
    return mainchain_all

if __name__ == '__main__':
    #def main():
    #argurements
        filename = sys.argv[1]
        #filename = ("file/%s" % pdb_name)
        cp = int(sys.argv[2])
        w = float(sys.argv[3])
        protein_name = sys.argv[1][-11:-7]
        start_time = time.time()
        print('reading data...')
        mainchain_all = Biopdb_read(filename)
        for num in range(len(mainchain_all)):
            res_num = mainchain_all[num][-1][-4] - mainchain_all[num][0][-4]+1
            start_resnum = mainchain_all[num][0][2]-1           
            mainchain = mainchain_all[num]        

            success_sum_array = np.zeros((res_num), dtype=np.int32)
            success_num_array = np.zeros((res_num, 40), dtype=np.int32)

            b = PDB_rot_read()
            b.read_pdbfile("file/rotamer.pdb")
            b.get_information_rot(res_names)
            b.mainchain_group(res_num, mainchain, res_names, start_resnum)
            rot_tuple = b.rot_tuple
            lines_rot = len(rot_tuple)
            resrotnum = b.resrotnum
            resatmnum = b.resatmnum
            res_atomnum_dict = b.res_atomnum_dict
            res_rotnum_dict = b.res_rotnum_dict
            rot_res = b.rot_res
            main_res_group = b.main_res_group
            main_res_name = b.main_res_name
            main_res_num = np.array(b.main_res_num, dtype = np.int32)
            main_res_atomnum = np.array(b.main_res_atomnum, dtype = np.int32)
            main_res_rotnum = np.array(b.main_res_rotnum, dtype = np.int32)
            
            c = Fit_it(
                res_num,
                main_res_name,
                main_res_group,
                rot_res,
                resrotnum,
                resatmnum,
                res_names,
                res_rotnum_dict,
                res_atomnum_dict,
                vdw,
                w)
            c.res_fit()
            res_rot_fited = c.res_rot_fited
            c.get_rc(res_atomnum_dict)
            c.get_rm(res_atomnum_dict)
            mrc = c.mrc
            mr = c.mr
            mrc_axyz = c.mrc_axyz
            mrc_d = c.mrc_d
            mr_d = c.mr_d
            atom_rc = c.atom_rc
            atom_r = c.atom_r
            atomnumc = c.atomnumc
            mainatom_xyz = reshape(mr_d, (4 * res_num, 3))
            mainatom_r = reshape(atom_r, (4 * res_num))
            print("time taken:" + ("%.2f" % float(time.time() - start_time)) + "seconds")
            print('caldist...')
            result = mc_rot(
                w,
                cp,
                res_num,
                atomnumc,
                main_res_atomnum,
                main_res_rotnum,
                mrc_axyz,
                mrc_d,
                mr_d,
                atom_r,
                atom_rc,
                mr,
                mrc,
                mainatom_r,
                mainatom_xyz,
                main_res_name,
                main_res_num,
                success_sum_array,
                success_num_array)
            cost_time = float(time.time() - start_time)
            file_name = 'result/result_20a'+str(w)+'0O.txt'
            with open (file_name,'a') as outfile:
                print("{}\t{:.2f}".format(res_num, result))
                outfile.write("{}\t{}\t{:.2f}\t{:.2f}s\n".format(protein_name,res_num, result, cost_time))
            print("time taken:" + ("%.2f" % float(time.time() - start_time)) + "seconds")
    #main()

