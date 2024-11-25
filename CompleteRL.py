# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 11:19:58 2024
author: Roger Chang 

- The stream of truth flows through its channels of mistakes.
"""

from tqdm import tqdm
import numpy as np

line_len = 60

"""
basic functions
"""

def splittingline(length):
    for i in range (length-1):
        print('-', end = '')
    print('-')
    return 0

def setdigit(x):
    return "{:10." + str(x) + "f}"

def loadFactors(path):
    F_raw = np.loadtxt(path, delimiter=",",
                        dtype=str, skiprows = 1)
    # replace all empty elems with 0
    F_raw[F_raw == ''] = 0
    # convert to float
    F = F_raw.astype(np.float32)
    return F

def factor2product(factors):
    # 1-D cases (only one factor)
    if(len(factors.shape) == 1):
        return factors
    else:
        product = np.array([1])
        for factor in range (factors.shape[0]):
            product = np.polymul(product, factors[factor, :])
        return product

def factors2TF(GN_factors, GD_factors, HN_factors, HD_factors):
    # forward gain: G(s) = GN(s)/GD(s)
    GN = factor2product(GN_factors)
    GD = factor2product(GD_factors)
    # feedback gain: H(s) = HN(s)/HD(s)
    HN = factor2product(HN_factors)
    HD = factor2product(HD_factors)
    # loop gain L = G * H
    LN = np.polymul(GN, HN)
    LD = np.polymul(GD, HD)
    return LN, LD

# frequency response: substitute s with jw
# return: real and imag part of freq resp function (func of w, discard j)
def freqResp(poly):
    poly_len = poly.shape[0]
    real_filter = np.zeros(poly_len)
    imag_filter = np.zeros(poly_len)
    for i in range (poly_len):
        deg = (poly_len-1) - i
        phase = deg % 4
        if(phase == 0):
            real_filter[i] =  1 
        elif(phase == 1):
            imag_filter[i] =  1
        elif(phase == 2):
            real_filter[i] = -1
        elif(phase == 3):
            imag_filter[i] = -1
        else:
            print("Unknown phase found.")
    real = poly * real_filter
    imag = poly * imag_filter

    return real, imag

"""
class: TF (transfer func)
"""
class TF:
    # initialize with coefficient csv paths
    def __init__(self, GN_path, GD_path, HN_path, HD_path):
        GN_factors = loadFactors(GN_path)
        GD_factors = loadFactors(GD_path)
        HN_factors = loadFactors(HN_path)
        HD_factors = loadFactors(HD_path)
        self.LN, self.LD = factors2TF(GN_factors, GD_factors, HN_factors, HD_factors)
    
    def OLPoles(self, precision):
        OLPs = np.roots(self.LD).astype(complex)
        print("Open loop poles")
        for pole in range (OLPs.shape[0]):
            print("  OLP. " + str(pole + 1) + ":  s = " + str(setdigit(precision).format(OLPs[pole])))
        splittingline(line_len)
        real = np.real(OLPs)
        imag = np.imag(OLPs)
        return real, imag
        
    def OLZeros(self, precision):
        OLZs = np.roots(self.LN).astype(complex)
        print("Open loop zeros")
        for zero in range (OLZs.shape[0]):
            print("  OLZ. " + str(zero + 1) + ":  s = " + str(setdigit(precision).format(OLZs[zero])))
        splittingline(line_len)
        real = np.real(OLZs)
        imag = np.imag(OLZs)
        return real, imag

    """
    asymptotes
    """
    # intersections
    def findInters(self, precision):
        OLZs = np.roots(self.LN)
        OLPs = np.roots(self.LD)
        m = OLZs.shape[0]
        n = OLPs.shape[0]
        asym_nums = n - m
        if(asym_nums > 0):
            sumPs = np.sum(OLPs)
            sumZs = np.sum(OLZs)
    
            intersect = (sumPs - sumZs)/asym_nums
            print("Intersection of asymptotes:")
            print("  Inters:  s = " + str(setdigit(precision).format(intersect)))
            splittingline(line_len)
            return [np.real(intersect)], [np.imag(intersect)]
        else:
            print("Intersection of asymptotes:     None")
            splittingline(line_len)
            return [], []
    # angles
    def findAngles(self, precision):
        OLZs = np.roots(self.LN)
        OLPs = np.roots(self.LD)
        n    = OLPs.shape[0]
        m    = OLZs.shape[0]
        asym_nums = n - m
        if(asym_nums > 0):
            k = np.arange(asym_nums)
            RL_angs_unit_pi = (2*k + 1)/asym_nums
            CR_angs_unit_pi =  2*k     /asym_nums
            print("Angles of asymptotes (k>0):")
            for i in range (asym_nums):
                print("  θr_" + str(i+1) + " = " 
                      + str(setdigit(precision).format(RL_angs_unit_pi[i])) 
                      + " π")
            print()
            print("Angles of asymptotes (k<0):")
            for j in range (asym_nums):
                print("  θc_" + str(j+1) + " = " 
                      + str(setdigit(precision).format(CR_angs_unit_pi[j]))
                      + " π")
            splittingline(line_len)
            return RL_angs_unit_pi * np.pi, CR_angs_unit_pi * np.pi
        else: 
            print("Angles of asymptotes (k>0): None")
            print("Angles of asymptotes (k<0): None")
            splittingline(line_len)
            return [], []
    
    """
    candidates of breakaway points
    """
    def cand_bwps(self, precision):
        # bwps occurrs at 1 + k * LN/LD = 0 has 2 or more identical roots
        # d/ds(1 + k * LN/LD) must equal to 0 when s = identical root value
        # d/ds(1 + k * LN/LD) = d/ds(LN/LD) = (LN' * LD - LD' * LN)/LD**2
        # LD ** 2 doesn't affect whether d/ds(LN/LD) = 0
        # solve (LN' * LD - LD' * LN) = 0
        dLN = np.polyder(self.LN, 1)
        dLD = np.polyder(self.LD, 1)
        eqLHS = np.polysub(np.polymul(dLN, self.LD), 
                           np.polymul(dLD, self.LN)) 
        cands = np.roots(eqLHS)
        if(cands.shape[0] > 0):
            print("Candidates of breakaway points:")
            for cand in range (cands.shape[0]):
                print("  cand. " + str(cand + 1) + ": s = " + str(setdigit(precision).format(cands[cand])))
            splittingline(line_len)
            return np.real(cands), np.imag(cands)
        else:
            print("Candidates of breakaway points: None")
            splittingline(line_len)
            return [], []
    
    """
    root locus
    """
    # forward_gain choosing range:
    # center: the gain that the highest degree term in k * LN
    #         has equal coefficeint in LD
    
    # define as k0 = LD[0]/LN[0]
    # sampling range: k0 * 10^k, k within [lLim, hLim]
    # with given sample amount
    
    # Root Locus - sampling roots in K>0
    def RL_poles(self, k, precision):
        k0 = self.LD[0]/self.LN[0]
        RL_CLPs  = np.zeros(1)
        print("Sampling root locus (k>0)...")
        for samp in tqdm(range (k.shape[0])):
            forward_gain = np.power(10, k[samp]) * k0
            GCL_CE = np.polyadd(forward_gain * self.LN, self.LD)
            samp_poles = np.roots(GCL_CE)
            # appending poles
            if(samp == 0):
                RL_CLPs = samp_poles
            else:    
                RL_CLPs = np.concatenate([RL_CLPs, samp_poles],
                                         axis=None)
        splittingline(line_len)
        RL_CLPs_real = np.real(RL_CLPs)
        RL_CLPs_imag = np.imag(RL_CLPs)
        return RL_CLPs_real, RL_CLPs_imag
        
    # Complementary Root Locus - sampling roots in K>0
    def CR_poles(self, k, precision):
        k0 = self.LD[0]/self.LN[0]
        CR_CLPs = np.zeros(1)
        print("Sampling complementary root locus (k<0)...")
        for samp in tqdm(range (k.shape[0])):
            forward_gain = -1 * np.power(10, k[samp]) * k0
            GCL_CE = np.polyadd(forward_gain * self.LN, self.LD)
            samp_poles = np.roots(GCL_CE)
            # appending poles
            if(samp == 0):
                CR_CLPs = samp_poles
            else:    
                CR_CLPs = np.concatenate([CR_CLPs, samp_poles], axis=None)
        splittingline(line_len)
        CR_CLPs_real = np.real(CR_CLPs)
        CR_CLPs_imag = np.imag(CR_CLPs)
        return CR_CLPs_real, CR_CLPs_imag

    """
    root locus - marginally stable points
    """
    def solveMarg(self, precision):
        """
        poly_real(w) = Re(poly(s=jw))
        poly_imag(w) = Im(poly(s=jw))
                
        eq1: P_real + k * Z_real = 0
        eq2: P_imag + k * Z_imag = 0
        
        <=> k = -P_real/Z_real = - P_imag/Z_imag
        <=> P_real / Z_real = P_imag / Z_imag
        <=> P_real * Z_imag = P_imag * Z_real
        """
        Z_real, Z_imag = freqResp(self.LN)
        P_real, P_imag = freqResp(self.LD)
        lhs = np.polymul(P_real, Z_imag)
        rhs = np.polymul(P_imag, Z_real)
        eq  = np.polysub(lhs, rhs)
        
        # only keep real w's
        marg_w = np.roots(eq)
        real_w = np.isreal(marg_w)
        marg_w = marg_w[real_w]
        marg_w = np.real(marg_w)
        # case: no real w's
        if(marg_w == []):
            print("Marginally-stable poles: None")
            return []
        else:            
            num_marg = int(marg_w.shape[0])
            print("Marginally-stable poles: " + str(num_marg) + " found.")
            marg_k = -1 * (np.polyval(P_real, marg_w) 
                           / np.polyval(Z_real, marg_w))

            for marg in range (num_marg):
                w = str(setdigit(precision).format(marg_w[marg]))
                k = str(setdigit(precision).format(marg_k[marg]))
                print("  pole " + str(marg + 1) + ":  s = " + w + "j" + ", k = " + k)
        splittingline(line_len)
        return marg_w, marg_k

        
"""
plotting boundary - xlim and ylim
"""

# xlim and ylim: maximum and minimum of open loop poles and zeros
# fixed aspect ratio: h:w = 3:4
# if xlim is set to -1, adjust automatically, minimum range is [-5, 5] for real axis
# er: expansion ratio
def setXYLims(xlim, er, 
              OLPs_real, OLPs_imag, 
              OLZs_real, OLZs_imag, 
              int_real,  int_imag, 
              cbwp_real, cbwp_imag, 
              msp_real,  msp_imag):   
    # aspect_ratio = h/w
    aspect_ratio = 0.73
    if(xlim == -1):
        points_real = np.concatenate([OLPs_real, OLZs_real, int_real, cbwp_real, msp_real])
        points_imag = np.concatenate([OLPs_imag, OLZs_imag, int_imag, cbwp_imag, msp_imag])
        real_max = np.max(np.absolute(points_real))
        imag_max = np.max(np.absolute(points_imag))

        if(imag_max > aspect_ratio * real_max):
            if(imag_max < aspect_ratio * 5 * er):
                imag_max = aspect_ratio * 5 * er
            imag_max = imag_max * er
            real_max = imag_max / aspect_ratio
        else:
            if(real_max < 5 * er):
                real_max = 5
            real_max = real_max * er
            imag_max = real_max * aspect_ratio
        return real_max, imag_max
    else:
        return xlim, xlim * aspect_ratio