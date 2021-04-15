# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 12:37:29 2018

@author: jm1e16
"""


def read_file(filename):
    """
    """

    data = [['Alpha', 'V', 'dE',
             'Xu', 'Xw', 'Zu', 'Zw', 'Zq', 'Mu', 'Mw', 'Mq',
             'Yv', 'Yp', 'Yr', 'Lv', 'Lp', 'Lr', 'Nv', 'Np', 'Nr',
             'Xde', 'Yde', 'Zde', 'Lde', 'Mde', 'Nde',
             'CXu', 'CXa', 'CZu', 'CLa', 'CLq', 'Cmu', 'Cma', 'Cmq',
             'CYb', 'CYp', 'CYr', 'Clb', 'Clp', 'Clr', 'Cnb', 'Cnp', 'Cnr',
             'CXde', 'CYde', 'CZde', 'Clde', 'Cmde', 'Cnde']]

    with open(filename, 'r') as f:
        lines = f.readlines()

    for i, l in enumerate(lines):
        if ('Rotating the elevator by' in l) or ('Rotating the flap by' in l):
            line = l
            de = line.split()[4][0:-2]
            line = lines[i+5]
        else: continue
        if 'Searching for zero-moment angle...' in line:
            alpha = line.split()[-1][6:-2]
            line = lines[i+8]
        else: continue
        if ('Calculating speed to balance the weight...' in line) and not('Found a negative lift' in line):
            V = line.split()[7]
            line = lines[i+26]
        else: continue

        if 'Xu' in line and 'Cxu' in line:
            Xu = line.split()[1]
            CXu = line.split()[3]
            line = lines[i+27]
        else: continue
        if 'Xw' in line and 'Cxa' in line:
            Xw = line.split()[1]
            CXa = line.split()[3]
            line = lines[i+28]
        else: continue
        if 'Zu' in line and 'Czu' in line:
            Zu = line.split()[1]
            CZu = line.split()[3]
            line = lines[i+29]
        else: continue
        if 'Zw' in line and 'CLa' in line:
            Zw = line.split()[1]
            CZa = line.split()[3]
            line = lines[i+30]
        else: continue
        if 'Zq' in line and 'CLq' in line:
            Zq = line.split()[1]
            CZq = line.split()[3]
            line = lines[i+31]
        else: continue
        if 'Mu' in line and 'Cmu' in line:
            Mu = line.split()[1]
            Cmu = line.split()[3]
            line = lines[i+32]
        else: continue
        if 'Mw' in line and 'Cma' in line:
            Mw = line.split()[1]
            Cma = line.split()[3]
            line = lines[i+33]
        else: continue
        if 'Mq' in line and 'Cmq' in line:
            Mq = line.split()[1]
            Cmq = line.split()[3]
            line = lines[i+38]
        else: continue

        if 'Yv' in line and 'CYb' in line:
            Yv = line.split()[1]
            CYb = line.split()[3]
            line = lines[i+39]
        else: continue
        if 'Yp' in line and 'CYp' in line:
            Yp = line.split()[1]
            CYp = line.split()[3]
            line = lines[i+40]
        else: continue
        if 'Yr' in line and 'CYr' in line:
            Yr = line.split()[1]
            CYr = line.split()[3]
            line = lines[i+41]
        else: continue
        if 'Lv' in line and 'Clb' in line:
            Lv = line.split()[1]
            Clb = line.split()[3]
            line = lines[i+42]
        else: continue
        if 'Lp' in line and 'Clp' in line:
            Lp = line.split()[1]
            Clp = line.split()[3]
            line = lines[i+43]
        else: continue
        if 'Lr' in line and 'Clr' in line:
            Lr = line.split()[1]
            Clr = line.split()[3]
            line = lines[i+44]
        else: continue
        if 'Nv' in line and 'Cnb' in line:
            Nv = line.split()[1]
            Cnb = line.split()[3]
            line = lines[i+45]
        else: continue
        if 'Np' in line and 'Cnp' in line:
            Np = line.split()[1]
            Cnp = line.split()[3]
            line = lines[i+46]
        else: continue
        if 'Nr' in line and 'Cnr' in line:
            Nr = line.split()[1]
            Cnr = line.split()[3]
            line = lines[i+49]
        else: continue

        if 'Xde' in line and 'CXde' in line:
            Xde = line.split()[1]
            CXde = line.split()[3]
            line = lines[i+50]
        else: continue
        if 'Yde' in line and 'CYde' in line:
            Yde = line.split()[1]
            CYde = line.split()[3]
            line = lines[i+51]
        else: continue
        if 'Zde' in line and 'CZde' in line:
            Zde = line.split()[1]
            CZde = line.split()[3]
            line = lines[i+52]
        else: continue
        if 'Lde' in line and 'CLde' in line:
            Lde = line.split()[1]
            Clde = line.split()[3]
            line = lines[i+53]
        else: continue
        if 'Mde' in line and 'CMde' in line:
            Mde = line.split()[1]
            Cmde = line.split()[3]
            line = lines[i+54]
        else: continue
        if 'Nde' in line and 'CNde' in line:
            Nde = line.split()[1]
            Cnde = line.split()[3]
            line = lines[i+55]
        else: continue


        data.append([alpha, V, de,
                     Xu, Xw, Zu, Zw, Zq, Mu, Mw, Mq,
                     Yv, Yp, Yr, Lv, Lp, Lr, Nv, Np, Nr,
                     Xde, Yde, Zde, Lde, Mde, Nde,
                     CXu, CXa, CZu, CZa, CZq, Cmu, Cma, Cmq,
                     CYb, CYp, CYr, Clb, Clp, Clr, Cnb, Cnp, Cnr,
                     CXde, CYde, CZde, Clde, Cmde, Cnde])

    return data


def file_write(filename):
    """
    """
    with open(filename, 'w') as file:
        for line in data:
            for s in line:
                file.write("{}, ".format(s))
            file.write("\n")


if __name__ == "__main__":
    data = read_file('T7-stability0.txt')
    print(data)
    file_write('T7-raw-out0.txt')

    data = read_file('T7-stability-dE.txt')
    print(data)
    file_write('T7-raw-out-dE.txt')

    data = read_file('T7-stability-dR.txt')
    print(data)
    file_write('T7-raw-out-dR.txt')


