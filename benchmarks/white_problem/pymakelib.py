# -*- coding: utf-8 -*-
'''
File: makelib.py
Project: MINICHEM
File Created: Tuesday, 2nd April 2019 11:24:00 am
Author: John Arul & Parth (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Last Modified: Friday, 5th April 2019 10:12:03 am
Modified By: John Arul & Parth Patel (arul@igcar.gov.in, parthpdpu@gmail.com)
-----
Copyright: IGCAR - 2019
'''
from scipy.io import FortranFile
import numpy as np
import fortranformat as ff
import re
'''
Program notes:
- It is confirmed that no gas has one temperature range (all have atleast two
  temperature ranges)
- Hence, while extrapolating the gases temperature coefficients, the program
  will write all three temperature ranges and corresponding coefficients
  including the extrapolated one.
- The condensed species temperature coeffcients are not extrapolted as the
  fortran program does not include it. In future these coefficients will be
  extrapolated.
- Formats used in this program to read fortran file:
  99001 FORMAT (4F10.3,a10)
  99002 FORMAT (a15,a65)
  99003 FORMAT (i2,1x,a6,1x,5(a2,f6.2),i2,f13.5,f15.3)
  99004 FORMAT (' ',a15,2x,a6,e15.6,2x,a65)
  99005 FORMAT (2F11.3,i1,8F5.1,2x,f15.3)
  99006 FORMAT (5D16.8/2D16.8,16x,2D16.8)
'''

f = open('thermo_python.inp', 'r')
op = open('thermo_chemical_database.txt', 'w')
tg_int = np.zeros([4])
sp_dict = {}
thermo_dict = {}


def filewrite_condensed_sp(sp_name, ntl, date, sp_dict,
                           phase, tl, mwt, thermo):
    op.write('{:20s}\t\t{:3.1f}\t\t{:10s}'.format(sp_name, ntl, date))
    for i12 in sp_dict[sp_name.replace('*', '')]:
        if i12[0] == '':
            i12[0] = '*'
        op.write('\t{:10s}\t{:3.1f}'.format(i12[0], i12[1]))
    # print(phase)
    op.write('{:5d}{:20.1f}{:20.1f}{:20.1e}'.format(phase, tl[0], tl[1], mwt))
    for j13 in range(thermo.shape[1]):
        for i13 in range(thermo.shape[0]):
            op.write('{:25.10e}'.format(thermo[i13, j13]))
    op.write('\n')


def filewrite_gas_sp(sp_name, ntl, date, sp_dict, phase,
                     total_temp_range, mwt, thermo):
    total_interval = ntl
    if total_interval <= 3:
        ntl = 3    # as all items are extrapolated
    op.write('{:20s}\t\t{:3.1f}\t\t{:10s}'.format(sp_name, ntl, date))
    for i12 in sp_dict[sp_name.replace('*', '')]:
        if i12[0] == '':
            i12[0] = '*'
        op.write('\t{:5s}\t\t{:3.1f}'.format(i12[0], i12[1]))
    # print(phase)
    op.write('{:5d}'.format(phase))
    for i23 in range(total_interval):
        op.write('{:15.1f}{:15.1f}'
                 .format(total_temp_range[sp_name.replace('*', '')][i23][0],
                         total_temp_range[sp_name.replace('*', '')][i23][1]))
    if total_interval < 3:
        op.write('{:25.1f}\t\t{:25.1f}\t'
                 .format(total_temp_range[sp_name.replace('*', '')][-1][1],
                         20000))
    op.write('{:20.8f}'.format(mwt))
    for j13 in range(thermo.shape[1]):
        for i13 in range(thermo.shape[0]):
            op.write('{:25.10e}'.format(thermo[i13, j13]))
    op.write('\n')

# skip to tables


i = 0
while(True):
    txt = f.readline()
    if re.match('thermo', txt):
        break

# read   header
ngl = 0
ns = 0
nall = 0
ifzm1 = 0
inew = 0
tinf = 1e6

txt = f.readline()
reader = ff.FortranRecordReader('4F10.3,a10')
arr = reader.read(txt)

tg_int = arr[:-1]
thedate = arr[-1]


def define_fill_thermo():
    fill = np.zeros(3)
    fill[:] = True
    thermo = np.zeros([9, 3])
    return fill, thermo

# sp_name = '   '


tl = np.zeros([2])

total_temp_range = {}
while True:

    fill, thermo = define_fill_thermo()

    hform = 0

    txt = f.readline()
    reader = ff.FortranRecordReader('a15,a65')
    arr = reader.read(txt)

    notes = arr[1].rstrip()
    sp_name = arr[0].rstrip()

    if (sp_name[:3] == 'END') or (sp_name[:3] == 'end') or (sp_name == 'Air'):

        ns = nall
        if ('ROD' not in sp_name) and ('rod' not in sp_name):

            break

    txt = f.readline()
    reader = ff.FortranRecordReader('i2,1x,a6,1x,5(a2,f6.2),i2,f13.5,f15.3')
    if txt == 'Air               \
            Mole:N2 78.084,O2 20.9476,Ar .9365,CO2 .0319.Gordon,1982.Reac\n':
        ns = nall
        break
    arr = reader.read(txt)
    ntl, date = arr[0:2]
    date = date.replace(' ', '')

    lis = []
    for i in range(5):
        j = 2 + (i * 2)
        lis.append([arr[j].replace(' ', ''), arr[j+1]])
    sp_dict.update({sp_name: lis})

    # phase of species 0: gas species
    #                  1: condensed species
    #                  2: contiguous species
    phase = arr[12]

    mwt = arr[13]     # molecule weight
    hform = arr[14]   # heat of formation

    if ntl == 0:

        if ns == 0:
            break

        nall = nall + 1

        txt = f.readline()
        reader = ff.FortranRecordReader('2F11.3,i1,8F5.1,2x,f15.3')
        arr = reader.read(txt)

        thermo[0, 0] = hform
        filewrite(sp_name, ntl, date, sp_dict, phase, tl, mwt, thermo)

    elif sp_name == 'Air':

        lis = [['N', 1.56168], ['O', 0.419590], ['AR', 0.009365], ['C',
               0.000319]]
        sp_dict.update({sp_name: lis})

    elif sp_name == 'e-':

        mwt = 5.48579903e-4
    if sp_name == 'ALBr':
        print('No of temp_interval = ', ntl)
    for i in range(ntl):

        txt = f.readline()
        reader = ff.FortranRecordReader('2F11.3,i1,8F5.1,2x,f15.3')
        arr = reader.read(txt)

        tl, ncoef, expn, hh = arr[0:2], arr[2:3], arr[3:10], arr[11]

        if sp_name not in list(total_temp_range.keys()):
            total_temp_range[sp_name] = []

        total_temp_range[sp_name].append(tl)
        txt1 = f.readline()
        txt2 = f.readline()

        reader = ff.FortranRecordReader('5D16.8/2D16.8,16x,2D16.8')
        templ = reader.read(txt1 + txt2)

        if phase == 0 and i > 3:
            # print('GO TO 400')
            break

        # For gas phase with three
        # temperature intervals, make
        # sp_name = *sp_name
        if phase <= 0:
            if tl[1] > (tg_int[3]-0.01):
                phase = -1
                sp_namee = '*' + sp_name
                sp_name = sp_namee
            if tl[0] >= tg_int[i + 1]:
                # print('GO TO 200')
                # GO TO 200
                break
            int = i

            fill[i] = False
        else:
            int = 0
            if i > 0:
                for k12 in range(7):
                    thermo[k12, 0] = 0
        for l in range(ncoef[0]):
            for k in range(7):
                if expn[l] == (k-2):
                    thermo[k, int] = templ[l]
                    #
                    # if sp_name == 'ALBr':
                    # print(thermo[k, int], k-2, l)
                    # breaker = True
                    break
            # if breaker:
            #     breaker = False
            #     break

        thermo[7, int] = templ[7]
        thermo[8, int] = templ[8]

        if phase > 0:
            nall = nall + 1
            if phase > ifzm1:
                inew = inew + 1
            else:
                inew = i

            filewrite_condensed_sp(sp_name, ntl, date, sp_dict, phase, tl, mwt,
                                   thermo)

    ifzm1 = phase
    if phase <= 0:
        inew = 0
        nall = nall + 1
        if phase <= 0 and ns == 0:

            ngl = ngl + 1
            if fill[2]:
                # print('hi ho')
                atms = 0
                for i1 in range(5):
                    if (sp_dict[sp_name][i1][0] == ' ') or (sp_dict[sp_name]
                                                            [i1][0] == 'E'):
                        # print('GO to 210')
                        break
                    # else:
                    atms = atms + sp_dict[sp_name][i1][1]

                aa = 2.5
                if atms > 1.9:
                    aa = 4.5
                if atms > 2.1:
                    aa = 3 * atms - 1.75
                ttl = tl[1]
                tx = ttl - tinf
                cpfix = 0
                templ[7] = 0
                templ[8] = 0
                dlt = np.log(ttl)
                for k13 in range(6, -1, -1):
                    kk = k13 - 2
                    # print(k13, kk)
                    if kk == 0:
                        cpfix = cpfix + thermo[k13, 1]
                        templ[7] = templ[7] + thermo[k13, 1]
                        templ[8] = templ[8] + thermo[k13, 1] * dlt
                        # print(templ[7])
                    else:
                        tex = ttl**kk
                        cpfix = cpfix + thermo[k13, 1] * tex
                        templ[8] = templ[8] + thermo[k13, 1] * tex / kk
                        if kk == -1:
                            # print(sp_name)
                            templ[7] = templ[7] + thermo[k13, 1] * dlt/ttl
                        else:
                            templ[7] = (templ[7]
                                        + thermo[k13, 1] * tex / (kk + 1))

                templ[1] = (cpfix-aa)/tx
                thermo[3, 2] = templ[1]
                templ[0] = cpfix - ttl * templ[1]
                thermo[2, 2] = templ[0]
                thermo[7, 2] = thermo[7, 1] + ttl * (templ[7] - templ[0] - 0.5
                                                     * templ[1] * ttl)
                thermo[8, 2] = (-templ[0] * dlt + thermo[8, 1] + templ[8] -
                                templ[1] * ttl)
        # filewrite(sp_name, ntl, date, sp_dict, phase, tl, mwt, thermo)
        filewrite_gas_sp(sp_name, ntl, date, sp_dict, phase, total_temp_range,
                         mwt, thermo)
        # WRITE (13) name,ntl,date,(sym(j),fno(j),j=1,5),ifaz,tl,mwt,
        #  &                thermo
    thermo_dict[sp_name.replace('*', '')] = thermo


op.close()
f.close()
