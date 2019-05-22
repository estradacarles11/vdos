import argparse
import numpy as np
import time

class Vdos():
    def __init__(self):
        
        self.def_natoms = 0
        self.def_lattct = 0
        self.def_lattctunits = "Ang"
        self.def_atcoordform = 'NotScaledCartesianBohr'
        self.def_atcoordorigin = [0, 0, 0]
        self.def_lxmax = 0
        self.def_lymax = 0
        self.def_lzmax = 0
        self.def_displ = 0.04
        self.def_displunits = "Bohr"
        self.def_blscale = "pi/a"
        
    def main(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str, help="<system-label>, FDF input file or Force Constant Matrix (.FC) input file.")
        parser.add_argument("-s", "--steps", type=int, help="Number of steps for the histogram.", default = 100)
        parser.add_argument("-nh", "--nkh", type=int, help="Number of divisions along reciprocal vector b1's direction.", default = 100)
        parser.add_argument("-nk", "--nkk", type=int, help="Number of divisions along reciprocal vector b2's direction.", default = 100)
        parser.add_argument("-nl", "--nkl", type=int, help="Number of divisions along reciprocal vector b3's direction.", default = 100)
        parser.add_argument("-o", "--output", type=str, help="Specify output file name.")
        self.args = parser.parse_args()
        
        if self.args.input.endswith('.FC'):
            self.syslabel = self.args.input[:-3:]
        elif self.args.input.endswith('.fdf'):
            self.syslabel = self.args.input[:-4:]
        else:
            self.syslabel = self.args.input
        self.sysfdf = "%s.fdf" % self.syslabel
        self.sysFC = "%s.FC" % self.syslabel
        
        self.steps = self.args.steps
        
        self.nk = [self.args.nkh, self.args.nkk, self.args.nkl]
        
        self.rnkh = xrange(self.nk[0])
        self.rnkk = xrange(self.nk[1])
        self.rnkl = xrange(self.nk[2])

        if self.args.output == None:
            self.sysvdos = "%s.vdos" % self.syslabel
        elif "." in self.args.output:
            self.sysvdos = self.args.output
        else:
            self.sysvdos = "%s.vdos" % self.args.output
            
        self.readfdf()
        fcvalues = self.readFC()
        fcmatrix = self.gendynmatrix(fcvalues)
        eigenvalues = self.geteigenval(fcmatrix)
        histogram = self.genhistogram(eigenvalues)
        self.writeoutput(histogram)
        
    def readfdf(self):
        try:
            inputfdf = open(self.sysfdf, 'r')
        except:
            print "File %s not found." % self.sysfdf
            exit()

        print ""
        print "Reading %s" % self.sysfdf
        
        self.natoms = self.def_natoms
        self.lattct = self.def_lattct
        self.lattctunits = self.def_lattctunits
        self.atcoordform = self.def_atcoordform
        self.atcoordorigin = self.def_atcoordorigin
        self.lxmax = self.def_lxmax
        self.lymax = self.def_lymax
        self.lzmax = self.def_lzmax
        self.displ = self.def_displ
        self.displunits = self.def_displunits
        self.blscale = self.def_blscale
        self.cell = np.zeros((3, 3))
        for i in xrange(3):
            self.cell[i][i] = 1
        self.ucell = np.zeros((3, 3))
        self.rcell = np.zeros((3, 3))
        
        scount = 1
        sc = "SuperCell_%s" % scount
        
        readacas = False
        readaco = False
        readbl = False
        readlp = False
        readlv = False
        
        fdflines = inputfdf.readlines()
        
        for line in fdflines:
            line = line.split()
            if len(line) == 2:
                if "NumberOfAtoms" in line[0]:
                    self.natoms = int(line[1])
                    if self.natoms > 0:
                        self.xa = np.zeros((3, self.natoms))
                        self.atindex = np.zeros(self.natoms)
                        self.xmass = np.zeros(self.natoms)
                    else:
                        print "Incorrect number of atoms in FDF."
                        exit()
                elif "%block" in line[0] and "LatticeParameters" in line[1]:
                    readlp = True
                elif "%block" in line[0] and "LatticeVectors" in line[1]:
                    readlv = True
                    vector = 0
                elif "%endblock" in line[0] and "LatticeVectors" in line[1]:
                    readlv = False
                elif "AtomicCoordinatesFormat" in line[0]:
                    self.atcoordform = str(line[1])
                elif "%block" in line[0] and "AtomicCoordinatesOrigin" in line[1]:
                    readaco = True
                elif "%endblock" in line[0] and "AtomicCoordinatesOrigin" in line[1]:
                    readaco = False
                elif "%block" in line[0] and "AtomicCoordinatesAndAtomicSpecies" in line[1]:
                    readacas = True
                    atom = 0
                elif "%endblock" in line[0] and "AtomicCoordinatesAndAtomicSpecies" in line[1]:
                    readacas = False
                elif sc in line[0]:
                    if scount == 1:
                        self.lxmax = int(line[1])
                    elif scount == 2:
                        self.lymax = int(line[1])
                    elif scount == 3:
                        self.lzmax = int(line[1])
                    scount += 1
                    sc = "SuperCell_%s" % scount
            elif len(line) == 3:
                if "LatticeConstant" in line[0]:
                    self.lattctunits = str(line[2])
                    if self.lattctunits == "Bohr":
                        self.lattct = float(line[1])
                    elif self.lattctunits == "Ang":
                        self.lattct = float(line[1]) / 0.529177
                elif "AtomicDispl" in line[0]:
                    self.displ = float(line[1])
                    self.displunits = str(line[2])
                elif readaco:
                    for i in xrange(3):
                        self.atcoord[i] = float(line[i])
                elif readlv:
                    for i in xrange(3):
                        self.cell[i][vector] = float(line[i])
                        """
                        if self.blscale == "ReciprocalLatticeVectors":
                            self.ucell[i][vector] = self.cell[i][vector] * self.lattct
                        """
                    vector += 1
            elif readlp:
                self.alp = float(line[0])
                self.blp = float(line[1])
                self.clp = float(line[2])
                self.alplp = float(line[3])
                self.betlp = float(line[4])
                self.gamlp = float(line[5])
                if int(self.alp) == 1 and int(self.blp) == 1 and int(self.clp) == 1 and int(self.alplp) == 90  and int(self.betlp) == 90 and int(self.gamlp) == 90:
                    for i in xrange(3):
                        self.ucell[i][i] = self.lattct

                self.cell = self.gencelllp(self.alp, self.blp, self.clp, self.alplp, self.betlp, self.gamlp)
                readlp = False
            elif readacas:
                for i in xrange(3):
                    self.xa[i][atom] = float(line[i])
                self.atindex[atom] = int(line[3])
                self.xmass[atom] = float(line[4])
                atom += 1
                    
        if self.natoms == 0:
            print "You MUST specify NumberOfatoms in FDF."
            exit()
        if self.lattct == 0:
            print "You MUST specify LatticeConstant in FDF."
            exit()
                    
        self.rnatoms = xrange(self.natoms)
        self.r3natoms = xrange(3 * self.natoms)
        self.rii = ["x", "y", "z"]
        self.r3 = xrange(3)
        self.rsign = ["negative", "positive", "average"]
        self.rlx = xrange(2 * self.lxmax + 1)
        self.rly = xrange(2 * self.lymax + 1)
        self.rlz = xrange(2 * self.lzmax + 1)
        self.ncells = (2 * self.lxmax + 1) * (2 * self.lymax + 1) * (2 * self.lzmax + 1)
        self.nnat = self.natoms * (2 * self.lxmax + 1) * (2 * self.lymax + 1) * (2 * self.lzmax + 1)
        self.lmax = [self.lxmax, self.lymax, self.lzmax]
        
        
        for ia in self.rnatoms:
            for ix in self.r3:
                self.xa[ix][ia] = self.xa[ix][ia] + self.atcoordorigin[ix]
        
        if self.atcoordform == "NotScaledCartesianBohr":
            pass
        elif self.atcoordform == "NotScaledCartesianAng":
            for ia in self.rnatoms:
                for ix in self.r3:
                    self.xa[ix][ia] = self.xa[ix][ia] / 0.529177
        elif self.atcoordform == "ScaledCartesian":
            for ia in self.rnatoms:
                for ix in self.r3:
                    self.xa[ix][ia] = self.xa[ix][ia] * self.lattct
        elif self.atcoordform == "ScaledByLatticeVectors":
            for ia in self.rnatoms:
                for ix in self.r3:
                    xac[ix] = self.xa[ix][ia]
                    self.xa[ix][ia] = 0
                for ix in self.r3:
                    for i in self.r3:
                        self.xa[ix][ia] += self.cell[ix][i] * xac[i]
        else:
            print "You must use one of the following coordinate scaling options:"
            print "\t- NotScaledCartesianBohr"
            print "\t- NotScaledCartesianAng"
            print "\t- ScaledCartesian"
            print "\t- ScaledByLatticeVectors"
        
        
        for i in self.r3:
            for ix in self.r3:
                self.cell[ix][i] = self.lattct * self.cell[ix][i]
        self.rcell = self.genrcell(self.cell)
        self.volume = np.linalg.det(self.cell)
        
        self.kpoint = np.zeros((3, self.nk[0], self.nk[1], self.nk[2]))
                
        for i in self.r3:    #actual xyz
            for nkh in self.rnkh:                #index for k point 1
                for nkk in self.rnkk:            #index for k point 2
                    for nkl in self.rnkl:        #index for k point 3
                        nkj = [nkh, nkk, nkl]
                        for j in self.r3:  #b1 -> b2 -> b3
                            self.kpoint[i][nkh][nkk][nkl] += nkj[j] * self.rcell[i][j] / self.nk[j]
        
        self.scell = np.zeros((3, 3))
        for ix in self.r3:
            for i in self.r3:
                self.scell[ix][i] = (2 * self.lmax[i] + 1) * self.cell[ix][i]
        
        """
        self.r = [0, 0, 0]
        iatom = 0
        for lx in self.rlx:
            for ly in self.rly:
                for lz in self.rlz:
                    for i in self.r3:
                        self.r(i) = lx * self.cell[i][0] + ly * self.cell[i][1] + lz * self.cell[i][2]
                    for i in self.rnatoms:
                        iatom += 1
                        imassc[iatom] = self.atindex[i]
                        for ix in self.r3:
                            self.sxa[ix][iatom] = self.xa[ix][i] + r[ix]
        
        
        self.i1 = self.natoms * (4 * self.lxmax * self.lymax * self.lzmax + 2 * (self.lxmax * self.lymax + self.lxmax * self.lzmax + self.lymax * self.lzmax) + self.lxmax + self.lymax + self.lzmax) + 1
        self.i2 = self.i1 + self.natoms - 1
        """
                        
                    

        print ""
        print "Atoms:\t%s" % self.natoms
        print "Lattice constant:\t%s %s" % (self.lattct, self.lattctunits)
        print "Unit cell:"
        print self.cell
        print "Reciprocal cell:"
        print self.rcell
        print "Atomic coordinates:"
        for i in self.rnatoms:
            print "%s\t%s\t%s\t%s\t%s (Bohr)"% (self.atindex[i], self.xmass[i], self.xa[0][i], self.xa[1][i], self.xa[2][i])
        print "lxmax:\t%s" % self.lxmax
        print "lymax:\t%s" % self.lymax
        print "lzmax:\t%s" % self.lzmax
        print "Cells in supercell:\t%s" % self.ncells
        print "Atomic displacement:\t%s %s" % (self.displ, self.displunits)
        print "Reciprocal space partitions\t%s" % self.nk
        print ""
        
        inputfdf.close()

    
    def gencelllp(self, alp, blp, clp, alplp, betlp, gamlp):
        
        alplp = alplp.np.deg2rad()
        betlp = betlp.np.deg2rad()
        gamlp = gamlp.np.deg2rad()
        
        xxx = (np.cos(alplp) - np.cos(betlp) * np.cos(gamlp)) /np.sin(gamplp)
        
        cell[0][0] = alp
        cell[1][0] = 0
        cell[2][0] = 0
        cell[0][1] = blp * np.cos(gamlp)
        cell[1][1] = blp * np.sin(gamlp)
        cell[2][1] = 0
        cell[0][2] = clp * cos(betlp)
        cell[1][2] = celp * xxx
        cell[2][2] = clp * (np.sin(betlp) * np.sin(betlp) - xxx*xxx)*0.5
        
        return cell

    def genrcell(self, A):
        
        B = np.zeros((3, 3))
        
        B[0][0] = A[1][1] * A[2][2] - A[2][1] * A[1][2]
        B[1][0] = A[2][1] * A[0][2] - A[0][1] * A[2][2]
        B[2][0] = A[0][1] * A[1][2] - A[1][1] * A[0][2]
        B[0][1] = A[1][2] * A[2][0] - A[2][2] * A[1][0]
        B[1][1] = A[2][2] * A[0][0] - A[0][2] * A[2][0]
        B[2][1] = A[0][2] * A[1][0] - A[1][2] * A[0][0]
        B[0][2] = A[1][0] * A[2][1] - A[2][0] * A[1][1]
        B[1][2] = A[2][0] * A[0][1] - A[0][0] * A[2][1]
        B[2][2] = A[0][0] * A[1][1] - A[1][0] * A[0][1]
        
        if self.blscale == "ReciprocalLatticeVectors":
            c = 2 * np.pi
        elif self.blscale == "pi/a":
            c = 1
        else:
            print "ERROR: Invalid value for BandLinesScale"
            exit()
        
        for i in xrange(3):
            ci = c / (A[0][i] * B[0][i] + A[1][i] * B[1][i] + A[2][i] * B[2][i])
            for j in xrange(3):
                B[j][i] = B[j][i] / ci
        
        return B
    
    def readFC(self):
        
        try:        
            inputFC = open(self.sysFC, 'r')   
        except:
            print "File %s not found." % self.sysFC
            exit()
        
        print "Reading %s" % self.sysFC
        print ""
        
        FCmatrix = np.loadtxt(inputFC, skiprows=1)
        
        inputFC.close()
                
        natoms = int(self.natoms)
        rnatoms = self.rnatoms
        r3 = self.r3
        rlx = self.rlx
        rly = self.rly
        rlz = self.rlz
        lxmax = self.lxmax
        lymax = self.lymax
        lzmax = self.lzmax
        ncells = self.ncells
        
        # index: [j][ij][lx][ly][lz][i][ii]
        pp = np.zeros((natoms, 3, 2 * lxmax + 1, 2 * lymax + 1, 2 * lzmax + 1, natoms, 3))
        pn = np.zeros((natoms, 3, 2 * lxmax + 1, 2 * lymax + 1, 2 * lzmax + 1, natoms, 3))
        phi0 = np.zeros((natoms, 3, 2 * lxmax + 1, 2 * lymax + 1, 2 * lzmax + 1, natoms, 3))
       
        #sign = 0            displacement sign (0 = neg, 1 = pos)
        #ij = 1               ??? (1, 2, 3)
        #i = 1               atom index in unit cell
        #j = 1                ??? (1 .. self.natoms)
        #lx = - self.lxmax   cell x coordinate in supercell
        #ly = - self.lymax   cell y coordinate in supercell
        #lz = - self.lzmax   cell z coordinate in supercell
        matrixc = 0
        
        for j in rnatoms:
            for ij in r3:
                for sign in r3:
                    for lx in rlx:
                        for ly in rly:
                            for lz in rlz:
                                for i in rnatoms:
                                    for ii in r3:
                                        if sign == 0:
                                            value = FCmatrix[matrixc][ii]
                                            pp[j][ij][lx][ly][lz][i][ii] = value
                                        elif sign == 1:
                                            value = FCmatrix[matrixc][ii]
                                            pn[j][ij][lx][ly][lz][i][ii] = value
                                        elif sign == 2:
                                            vpp = pp[j][ij][lx][ly][lz][i][ii]
                                            vpn = pn[j][ij][lx][ly][lz][i][ii]
                                            phi0[j][ij][lx][ly][lz][i][ii] = (vpp + vpn) / 2
                                            if ii == 0:
                                                matrixc -= 1
                                        #print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (matrixc, value, j, ij, rsign[sign], lx, ly, lz, i, rii[ii])
                                    matrixc += 1
        
        
        #phibar(1)  (it should be hermitian)
        
        # index: [j][ij][lx][ly][lz][i][ii]
        phi = np.zeros((natoms, 3, 2 * lxmax + 1, 2 * lymax + 1, 2 * lzmax + 1, natoms, 3))
        
        for i in rnatoms:
            for j in rnatoms:
                for ii in r3:
                    for ij in r3:
                        for lx in rlx:
                            for ly in rly:
                                for lz in rlz:
                                    vphip = phi0[j][ij][lx][ly][lz][i][ii]
                                    vphin = phi0[i][ii][lxmax - 1 - lx][lymax - 1 - ly][lzmax - 1 - lz][j][ij]
                                    phi[j][ij][lx][ly][lz][i][ii] = (vphip + vphin) / 2
        
        """
        goto 100
        """
        
        # index [ii][ij][j]
        zero = np.zeros((3, 3, natoms))
        
        for j in rnatoms:
            for ii in r3:
                for ij in r3:
                    zero[ii][ij][j] = 0
                    for lx in rlx:
                        for ly in rly:
                            for lz in rlz:
                                zero[ii][ij][j] = zero[ii][ij][j] + phi[j][ij][lx][ly][lz][i][ii]
                    zero[ii][ij][j] = zero[ii][ij][j] / (natoms * ncells)
                    
        # index [ii][ij]
        zeroo = np.zeros((3, 3))
        
        for ii in r3:
            for ij in r3:
                zeroo[ii][ij] = 0
                for j in rnatoms:
                    zeroo[ii][ij] = zeroo[ii][ij] + zero[ii][ij][j]
                zeroo[ii][ij] = zeroo[ii][ij] / natoms
                
        #phibar(2) (it IS hermitian)
        
        # index: [j][ij][lx][ly][lz][i][ii]
        for i in rnatoms:
            for j in rnatoms:
                for ii in r3:
                    for ij in r3:
                        correct = (zeroo[ii][ij] + zeroo[ij][ii]) /2 - (zero[ii][ij][j] + zero[ij][ii][i])
                        for lx in rlx:
                            for ly in rly:
                                for lz in rlz:
                                    phi[j][ij][lx][ly][lz][i][ii] = phi[j][ij][lx][ly][lz][i][ii] + correct
                                    
        return phi
        
        
    def gendynmatrix(self, phi):
        
        print "Generating matrix"
                
        nk = self.nk
        rnkh = self.rnkh
        rnkk = self.rnkk
        rnkl = self.rnkl
        r3 = self.r3
        kpoint = self.kpoint
        rnatoms = self.rnatoms
        rlx = self.rlx
        rly = self.rly
        rlz = self.rlz
        xa = self.xa
        cell = self.cell
        scell = self.scell
        r3natoms = self.r3natoms
        xmass = self.xmass
        natoms = int(self.natoms)
                
        ### vibrator.f line 452+

        dc3 = natoms * 3
        
        qz3 = np.zeros(3)
        rz3 = np.zeros(3)
        r1z3 = np.zeros(3)
        dcz = np.zeros((dc3, dc3, 2))
        dd = np.zeros((nk[0], nk[1], nk[2], dc3, dc3))
        

        
        for ikx in rnkh:                   #for ikx in rnkh:      for ikx in [0]:
            for iky in rnkk:               #for iky in rnkk:      for iky in [0]:
                for ikz in rnkl:           #for ikz in rnkl:      for ikz in [0]:
                    q = qz3[:]
                    for x in r3:
                        q[x] = kpoint[x][ikx][iky][ikz]
                    dc = dcz[:]
                    #dd = ddz[:]
                    for i in rnatoms:
                        for j in rnatoms:
                            for lx in rlx:
                                    for ly in rly:
                                        for lz in rlz:
                                            r = rz3[:]
                                            for jj in r3:
                                                r[jj] = xa[jj][i] - xa[jj][j] + lx * cell[jj][0] + ly * cell[jj][1] + lz * cell[jj][2]
                                            dmin = 1.0e+08
                                            for llx in r3:
                                                for lly in r3:
                                                    for llz in r3:
                                                        r2 = 0
                                                        r1 = r1z3[:]
                                                        for jj in r3:
                                                            r1[jj] = llx * scell[jj][0] + lly * scell[jj][1] + llz * scell[jj][2]
                                                            r2 = r2 + (r1[jj] + r[jj])**2
                                                        absrd = abs(r2 - dmin)
                                                        if absrd > 0.0001:
                                                            if r2 < dmin:
                                                                neq = 1
                                                                dmin = r2
                                                                qr = [q[0] * (r1[0] + r[0]) + q[1] * (r1[1] + r[1]) + q[2] * (r1[2] + r[2])]
                                                        elif absrd < 0.0001:
                                                            neq += 1
                                                            qr.append(q[0] * (r1[0] + r[0]) + q[1] * (r1[1] + r[1]) + q[2] * (r1[2] + r[2]))
                                            for inn in xrange(neq):
                                                phase = [np.cos(qr[inn]), np.sin(qr[inn])]
                                                for ii in r3:
                                                    for ij in r3:
                                                        ix = (i - 1) * 3 + ii
                                                        jx = (j - 1) * 3 + ij
                                                        dcsum = [0, 0]
                                                        index = 0
                                                        for ph in phase:
                                                            dcsum[index] = phi[j][ij][lx][ly][lz][i][ii] * ph / neq
                                                            index += 1
                                                        for i in [0, 1]:
                                                            dc[ix][jx][i] += dcsum[i]
                    
                    #print dc
                    
                    for ix in r3natoms:
                        for jx in r3natoms:
                            modix = np.sqrt((ix - 1)**2 + 9)
                            modjx = np.sqrt((jx - 1)**2 + 9)
                            i = int((ix - modix) / 3 + 1)
                            j = int((jx - modjx) / 3 + 1)
                            for k in [0, 1]:
                                dc[ix][jx][k] = dc[ix][jx][k] / np.sqrt(xmass[i] * xmass[j])
                                
                    for ix in r3natoms:
                        for jx in r3natoms:
                            dd[ikx][iky][ikz][ix][jx] = dc[jx][ix][1]  #np.imag(dc[jx][ix])
                            dd[ikx][iky][ikz][jx][ix] = dc[jx][ix][0]   #np.real(dc[jx][ix])
                    
                    #print ikx, iky, ikz
                    #print q
                    #print ""
                    #print dd[ikx][iky][ikz]
                    #print ""
                    
        return dd
    
    
    def geteigenval(self, dd):
        
        natoms = int(self.natoms)
        r3natoms = self.r3natoms
        nk = self.nk
        rnkh = self.rnkh
        rnkk = self.rnkk
        rnkl = self.rnkl
        kpoint = self.kpoint
        
        eigenvalues = np.zeros((nk[0] * nk[1] * nk[2] * 3 * natoms, 2))
        
        self.minen = 0
        self.maxen = 0
        
        count = 0
        for ikx in rnkh:                   #for ikx in rnkh:      for ikx in [0]:
            for iky in rnkk:               #for iky in rnkk:      for iky in [0]:
                for ikz in rnkl:           #for ikz in rnkl:      for ikz in [0]:
                    r = np.sqrt(kpoint[0][ikx][iky][ikz]**2 + kpoint[1][ikx][iky][ikz]**2 + kpoint[2][ikx][iky][ikz]**2)
                    for i in r3natoms:
                        eigenvalues[count][0] = r
                        eigenvalues[count][1] = np.linalg.eigvalsh(dd[ikx][iky][ikz])[i]
                        if eigenvalues[count][1] < self.minen:
                            self.minen = eigenvalues[count][1]
                        elif eigenvalues[count][1] > self.maxen:
                            self.maxen = eigenvalues[count][1]
                        count += 1
        print ""
        print eigenvalues
        
        return eigenvalues
    
    def genhistogram(self, eigenvalues):
        
        natoms = int(self.natoms)
        nk = self.nk
        
        maxenergy = self.maxen / 8065.73
        minenergy = self.minen / 8065.73
        
        k = np.zeros(nk[0] * nk[1] * nk[2] * 3 * natoms)
        energy = np.zeros(nk[0] * nk[1] * nk[2] * 3 * natoms)
        
        for i in xrange(nk[0] * nk[1] * nk[2] * 3 * natoms):
            k[i] = eigenvalues[i][0]
            energy[i] = eigenvalues[i][1] / 8065.73
        
        step = (maxenergy - minenergy) / self.steps
        binss = [minenergy]
        while binss[-1] < maxenergy:
            binss.append(binss[-1] + step)
        
        hist, bin_edges = np.histogram(energy, bins = binss)
        
        l = range(len(hist))
        histogram = np.zeros((len(hist),2))
        for i in l:
            histogram[i][0] = bin_edges[i]
            histogram[i][1] = hist[i] / step / (nk[0] * nk[1] * nk[2])
        
        return histogram
    
    def writeoutput(self, histogram):
        output = open(self.sysvdos, 'w')
        for i in histogram:
            output.write('%s\t%s\n' % (i[0], i[1]))
        output.close()
            


starttime = time.asctime(time.localtime(time.time()))
vdos = Vdos()        
vdos.main()

endtime = time.asctime(time.localtime(time.time()))
print "Start:", starttime
print "End:", endtime
