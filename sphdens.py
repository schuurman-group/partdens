#!/usr/bin/env python
"""
Script to generate radial density matrices for atoms

To obtain radial densities, the overlap can be expressed as angular
and radial quantities

    rho = Tr(DS) = sum_ij D_ij S_ij = sum_ij D_ij A_ij R_ij.

Integrating out angular terms leads to all terms within the same shell
having the same R_ij. Thus, the radial density matrix can be reduced
from n_basis x n_basis to n_shell x n_shell

    D_r,ij' = D_ij A_ij,
    D_r,ij  = T^T D_r,ij' T,

where T is an n_basis x n_shell matrix which sums over all basis
functions in the same shell.

The elements of A are

    A_ij = N_i N_j int_0^2pi dy cos^(ai+aj) y sin^(bi+bj) y
           int_0^pi dz sin^(1+ai+aj+bi+bj) z cos^(ci+cj) z,

where N_i and N_j satisfy A_ii = 1. The integrals yield 0 if
xi + xj is odd for x in {a, b, c}. Otherwise,

    A_ij = sqrt[(2li + 1)!!(2lj + 1)!!] / (li + lj + 1)!! *
           prod_x (xi + xj - 1)!! / sqrt[(2xi - 1)!!(2xj - 1)!!],

where li = sum_x xi.

The angular component is integrated out if the flag '-a' is included,
otherwise the density matrix is simply symmetrized and reduced.
"""
import sys
import numpy as np


def read_stub(string):
    """Reads the element label and filename from a file stub."""
    lbl, chgstr = string.split('_')
    chg = chgstr.count('p') - chgstr.count('m')
    return lbl, chg


def read_mos(fname, met='scf'):
    """Reads the MOs from a GAMESS data file."""
    if met == 'scf':
        met_string = 'ROHF'
    elif met == 'mp2':
        met_string = 'MP2 NATURAL ORBITALS'

    with open(fname, 'r') as f:
        dat = f.readlines()

    idat = 0
    for line in dat:
        if '$DATA' in line:
            idat += 4
            break
        else:
            idat += 1

    lbl = []
    while len(dat[idat].split()) > 0:
        lshell, ncart = dat[idat].split()
        lbl.append(lshell)
        idat += int(ncart) + 1

    vec = []
    raw_mos = ''
    last_i = -1
    read = False
    for line in dat:
        if '$VEC' in line:
            continue
        elif met_string in line:
            read = True
        elif read and '$END' in line:
            vec.append([float(raw_mos[j:j+15]) for j
                        in range(0, len(raw_mos), 15)])
            break
        elif read:
            i = int(line[:2]) - 1
            if i > last_i:
                if raw_mos != '':
                    vec.append([float(raw_mos[j:j+15]) for j
                                in range(0, len(raw_mos), 15)])
                raw_mos = line[5:].replace('\n', '')
                last_i = 1*i
            else:
                raw_mos += line[5:].replace('\n', '')

    return lbl, np.array(vec)


def read_out(fname, met='scf'):
    """Reads the basis and occupations from a GAMESS output file."""
    with open(fname, 'r') as f:
        dat = f.readlines()

    for line in dat:
        if '$BASIS REQUESTS READING THE' in line:
            bas = line.split('"')[1].strip()
            break

    if met == 'scf':
        for line in dat:
            if 'NUMBER OF CARTESIAN' in line:
                occ = np.zeros(int(line.split()[-1]))
            elif 'OCCUPIED ORBITALS (ALPHA)' in line:
                occ[:int(line.split()[-1])] += 1
            elif 'OCCUPIED ORBITALS (BETA )' in line:
                occ[:int(line.split()[-1])] += 1
            elif 'TOTAL NUMBER OF ATOMS' in line:
                break
    if met == 'mp2':
        occ = []
        read = False
        for line in dat:
            if 'MP2 NATURAL ORBITAL OCCUPATION' in line:
                read = True
            elif read and 'PRINCIPAL MP2' in line:
                occ = np.array(occ, dtype=float)
                break
            elif read:
                occ += [float(line[j:j+15]) for j in range(1, len(line)-1, 15)]

    return bas, occ


def build_mats(labels, ang=True):
    """Returns the angular overlap matrix and the number of AOs per
    shell based on a list of shell labels."""
    nlbl = len(labels)

    # find T and powers of x, y, and z
    ncart = np.empty((0, 3), dtype=int)
    nblk = np.empty(0, dtype=int)
    nao = 0
    for i, l in enumerate(labels):
        n, eq = get_nxyz(l)
        ncart = np.vstack((ncart, n))
        if ang:
            nblk = np.hstack((nblk, len(n)))
        else:
            nblk = np.hstack((nblk, eq))

        nao += len(n)

    if ang:
        # fill in the lower triangular of A
        amat = np.zeros((nao, nao))
        for i in range(nao):
            amat[i,i] = 1
            for j in range(i):
                ni = ncart[i]
                nj = ncart[j]
                if np.all(ni == nj):
                    amat[i,j] = 1
                elif np.all((ni + nj) % 2 == 0):
                    li = np.sum(ni)
                    lj = np.sum(nj)
                    nij = dfac(ni + nj - 1) / np.sqrt(dfac(2*ni - 1)*dfac(2*nj - 1))
                    amat[i,j] = np.prod(nij)
                    if li != lj:
                        amat[i,j] *= np.sqrt(dfac(2*li + 1)*dfac(2*lj + 1))
                        amat[i,j] /= dfac(li + lj + 1)

        # fill in the upper triangular of A
        iu = np.triu_indices(nao, k=1)
        il = (iu[1], iu[0])
        amat[iu] = amat[il]
        return amat, nblk
    else:
        return 1, nblk


def get_nxyz(lbl):
    """Returns the gamess-ordered x, y and z powers for a shell.

    Given l = i + j + k, i >= j >= k, the larger number is assigned to
    x, y, z in that order. j and i are decreased until the equality is
    no longer satisfied. The ordering of g and h orbitals need to be
    double-checked for consistency with gamess output.
    """
    lxyz = dict(S=0, P=1, D=2, F=3, G=4, H=5, I=6)
    n = lxyz[lbl]
    if n == 0:
        return [[0, 0, 0]], [1]

    plist = [[n, 0, 0], [0, n, 0], [0, 0, n]]
    neq = [3]
    for i in range(n, 0, -1):
        for j in range(n - i, 0, -1):
            k = n - i - j
            if k > j:
                break
            elif j > i:
                continue
            elif i == j and j == k:
                plist += [[i, i, i]]
                neq += [1]
            elif i == j:
                plist += [[i, i, k], [i, k, i], [k, i, i]]
                neq += [3]
            elif j == k:
                plist += [[i, j, j], [j, i, j], [j, j, i]]
                neq += [3]
            else:
                plist += [[i, j, k], [i, k, j], [j, i, k],
                          [j, k, i], [k, i, j], [k, j, i]]
                neq += [6]

    return plist, neq


def dfac(x):
    """Returns the double factorial of an integer or array of integers."""
    xi = np.copy(x)
    xfac2 = np.ones_like(x)
    while np.any(xi > 0):
        xfac2[xi > 0] *= xi[xi > 0]
        xi -= 2

    return xfac2


def write_densmat(f, dens, lbl, bas, chg, neq=None, offd=None):
    """Writes the symmetrized density matrix to an output file."""
    nblk = len(dens)
    nnum = (nblk // 5)*[5] + [nblk % 5]

    # write atom, basis and charge
    if chg == 0:
        f.write('{:3s}{:8s}{:3d}\n'.format(lbl, bas, chg))
    else:
        f.write('{:3s}{:8s}{:+3d}\n'.format(lbl, bas, chg))

    if neq is not None:
        # write number of equivalent blocks
        ipr = 0
        for i, n in enumerate(nnum):
            fmt = '  {:3d}' + n*'{:15d}' + '\n'
            f.write(fmt.format(i+1, *neq[ipr:ipr+n]))
            ipr += n

    if offd is not None:
        # write off-diagonal elements
        ipr = 0
        for i, n in enumerate(nnum):
            fmt = '  {:3d}' + n*'{:15.8E}' + '\n'
            f.write(fmt.format(i+1, *offd[ipr:ipr+n]))
            ipr += n

    # write density matrix
    for i, d in enumerate(dens):
        ipr = 0
        for j, n in enumerate(nnum):
            fmt = '{:2d}{:3d}' + n*'{:15.8E}' + '\n'
            f.write(fmt.format((i+1) % 100, j+1, *d[ipr:ipr+n]))
            ipr += n


def main():
    """The main routine."""
    if '-mp2' in sys.argv:
        method = 'mp2'
        sys.argv.pop(sys.argv.index('-mp2'))
    else:
        method = 'scf'

    integ_ang = '-a' in sys.argv
    if integ_ang:
        sys.argv.pop(sys.argv.index('-a'))

    # read in gamess data
    stub = sys.argv[1].replace('.dat', '')
    atm, chg = read_stub(stub)
    l, c = read_mos(stub + '.dat', met=method)
    basis, n = read_out(stub + '.out', met=method)
    nao = len(c)

    # get A and numbers of equivalent AOs
    ang, neqiv = build_mats(l, ang=integ_ang)
    nblk = len(neqiv)
    tra = np.zeros((nblk, nao), dtype=int)
    iao = 0
    for i, eq in enumerate(neqiv):
        tra[i,iao:iao+eq] = 1
        iao += eq

    # form D and reduce it to D_r
    da = ang*np.dot(c.T*n, c)
    dr = np.dot(tra, np.dot(da, tra.T))

    if integ_ang:
        neqiv = None
        od = None
        sumdr = np.sum(dr)
    else:
        od = np.zeros(nblk)
        iao = 0
        for i, eq in enumerate(neqiv):
            dr[i,i] = np.sum(np.diag(da[iao:iao+eq,iao:iao+eq]))
            od[i] = np.sum(np.tril(da[iao:iao+eq,iao:iao+eq], k=-1))
            iao += eq

        sumdr = np.sum(dr) + 2*np.sum(od)

    # check values of D_r
    if not np.isclose(np.sum(da), sumdr):
        raise ValueError('Sum of D_r not conserved ' +
                         '({:10.3e},{:10.3e})'.format(np.sum(da), sumdr))
    if not integ_ang:
        if not np.isclose(np.trace(da), np.trace(dr)):
            raise ValueError('Trace of D_r not conserved ' +
                             '({:10.3e},{:10.3e})'.format(np.trace(da), np.trace(dr)))

    # output D_r
    write_densmat(sys.stdout, dr, atm, basis, chg, neq=neqiv, offd=od)


if __name__ == '__main__':
    main()
