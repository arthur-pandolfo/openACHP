# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 22:25:20 2015

The purpose of this file is to check the compatibility of the relations from
Patek and Klomfar, Int J Refrig, Vol 20, pp 566-578 (2006) with the CoolProp
functions for water, and compare to a known set of code using these relations.

Units used:
    P : float
        pressure, [bar]
    T : float
        temperature, [K]
    h : float
        enthalpy, [J/kg]
    s : float
        entropy, [J/kg-K]
    rho : float
        density, [kg/m3]
    MW : float
        molecular weight, [kg/mol]

Crystallization curve
---------------------
    libr_props.crystallization_data_T :
        temperature, degrees C
    libr_props.crystallization_data_x :
        mass fraction, kg/kg

Source: Boryta, D.A., 1970, "Solubility of Lithium Bromide in Water between
-50 and +100 C (40 to 70% LiBr)", J Chem Eng Data, Vol 15, No 1, pp 142-144.

@author: nfette
"""

import CoolProp as CP
from CoolProp.CoolProp import PropsSI
from CoolProp import AbstractState, constants
from hw2_1 import CelsiusToKelvin as C2K, KelvinToCelsius as K2C
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import numpy as np

MW_LiBr = 0.08685 # kg/mol
MW_H2O = 0.018015268 # kg/mol
xmax = 0.7
pwater = AbstractState("HEOS","water")
# The LiBr enthalpy is zero at 293.15 K, 1 atm, per
# http://www.coolprop.org/fluid_properties/Incompressibles.html#general-introduction
# We need to evaluate water enthalpy relative to that, but it is not built-in.
# http://www.coolprop.org/coolprop/HighLevelAPI.html#reference-states
T_ref = 20
P_ref = 101325

pwater.update(CP.PT_INPUTS,P_ref,C2K(T_ref))
h_w_ref = pwater.hmass()
h_w_molar_ref = pwater.hmolar()

def mole2massFraction(x):
    """input: mole fraction, x, of LiBr"""
    return x * MW_LiBr / (x * MW_LiBr + (1 - x) * MW_H2O)
    
def molefraction(w):
    """input: mass fraction, w, of LiBr"""
    return (w / MW_LiBr) / (w/MW_LiBr + (1 - w) / MW_H2O)

def thetaFun(T,x,Tderiv=False,Xderiv=False):
    a = [-2.41303e2, 1.91750e7, -1.75521e8, 3.25432e7,
         3.92571e2, -2.12626e3, 1.85127e8, 1.91216e3] # [K]
    m = [3,4,4,8,1,1,4,6]
    n = [0,5,6,3,0,2,6,0]
    t = [0,0,0,0,1,1,1,1]
    T_c = 647.096 # [K]
    TK = T
    x_N = molefraction(x)
    
    s = 0
    for i in range(8):
        s=s+a[i]* x_N**m[i] * (0.4-x_N)**n[i] * (TK/T_c)**t[i]
    Theta=TK-s
    result = (Theta,)
    
    # Next compute the derivative wrt T
    if Tderiv:
        s = 0
        for i in range(4,8):
            s = s + a[i] * x_N ** m[i] * (0.4 - x_N) ** n[i] / T_c
        dThdT = 1 - s
        result = result + (dThdT,)

    # Next compute the derivative wrt x
    if Xderiv:
        s1, s2 = 0, 0
        for i in range(8):
            s1 = s1 + a[i] * m[i] * x_N ** (m[i] - 1) * (0.4 - x_N) ** n[i] \
                 * (TK / T_c) ** t[i]
            if (n[i] > 0):
                s2 = s2 + a[i] * x_N ** m[i] * n[i] * (0.4 - x_N) ** (n[i] - 1) \
                     * (TK / T_c) ** t[i]
        dThdx = s2 - s1
        result = result + (dThdx,)
        
    return result

def pressure(T,x):
    """Return pressure above a lithium bromide-water mixture
    at the given temperature and composition using the formulation
    presented by
    Patek and Klomfar, Int. J. Refrig., Vol 29, pp 566-578 (2006)
    
    Notes: "above" the mixture: is completely water vapor. So there are only
    two relevant properties to find equilibrium vapor pressure (?).

    Units: T [K]
           x = mass fraction LiBr
           P [bar]
           
    Based on Table 4 and Equation (1) in reference.
    """
    
    Theta, = thetaFun(T,x)
    Q = 0.0
    #print("Theta, Q = {}, {}".format(Theta, Q))
    pwater.update(CP.QT_INPUTS, Q, Theta)
    pressurePa = pwater.p()
    #print("pressurePa = {}".format(pressurePa)) # ()
    pressureBar = pressurePa * 1e-5
    return pressureBar
    
def objective_T(T,*ThetaX):
    Theta,x= ThetaX
    ThetaOut,dThdT = thetaFun(T,x,Tderiv=True)
    #print("Theta,x,T,ThetaOut,dThdT = {}, {}, {}, {}, {}"
    #      .format(Theta,x,T,ThetaOut,dThdT))
    return [(ThetaOut - Theta) ** 2, 2 * (ThetaOut - Theta) * dThdT]
    
def temperature(P,x):    
    """T_LiBrH2O returns the temperature of a lithium bromide-water mixture at
    the given the pressure and composition using the formulation presented by
    Patek and Klomfar, Int. J. Refrig., Vol 29, pp 566-578 (2006)
    
    Notes: "above" the mixture: is completely water vapor. So there are only
    two relevant properties to find equilibrium vapor pressure (?).

    Units:
        T [K]
        x = mass fraction LiBr
        P [bar]
    """
    global soln
    # Just call a solver on the previously defined pressure function
    # Trick is to get the guess temperature on the high side.
    # Does this constraint do anything? Maybe not.
    #cons = ({"type": "ineq",
    #         "fun": lambda T: np.array(P - pressure(T,x)),
    #        },)
    #soln = minimize(objective_T, guess, constraints=cons, args=(P,x,))
    Q = 0.0
    pressurePa = P * 1e5
    pwater.update(CP.PQ_INPUTS, pressurePa, Q)
    theta = pwater.T() # K
    guess = (theta,)
    soln = minimize(objective_T, guess, args=(theta,x),jac=True, bounds=[(0.,647.)])
    #print("Success, message: {}, {}".format(soln.success, soln.message))
    return soln.x[0]

def objective_x(x,*TTheta):
    #print("T, Px = {}, {}".format(T,Px))
    T,Theta = TTheta
    ThetaOut,dThdx = thetaFun(T, x, Xderiv=True)
    f,fprime = (ThetaOut - Theta) ** 2, 2 * (ThetaOut - Theta) * dThdx
    #print("x,ThetaOut,dThdx,f,f' = {}, {}, {}, {}, {}"
    #      .format(x,ThetaOut,dThdx,f,fprime))
    return [f,fprime]
    
def massFraction(T,P,guess=(0.5,)):
    """Returns the composition of a lithium bromide-water mixture at
    the given the temprature and pressure using the formulation presented by
    Patek and Klomfar, Int. J. of Refrigeration, Vol 29, pp. 566-578, (2006)

    Notes: "above" the mixture: is completely water vapor. So there are only
    two relevant properties to find equilibrium vapor pressure (?).

    Args
    ----
        T [K]
            Temperature
        P [bar]
            Pressure
    
    Outputs
    -------
        x [kg/kg]
            Mass fraction LiBr
    """
    global soln2
    # Just call a solver on the previously defined pressure function
    # Trick required? Guess mass fraction not too high nor low.
    # Does this constraint do anything? Maybe not.
    #cons = ({"type": "ineq",
    #         "fun": lambda x: np.array(pressure(T,x) - P),
    #        },)
    #soln = minimize(objective_x, guess, constraints=cons, args=(T,P))
    pressurePa = P * 1e5
    Q = 0.0
    pwater.update(CP.PQ_INPUTS, pressurePa, Q)
    theta = pwater.T() # [K]
    #print("T,P,Theta = {}, {}, {}".format(T,P,theta))
    soln2 = minimize(objective_x, guess, args=(T,theta),jac=True,bounds=[(0,1)])
    #print("Success, message: {}, {}".format(soln2.success, soln2.message))
    return soln2.x[0]

def massSpecificEnthalpy(T,x):
    """Inputs:  T = Temperature / [Kelvin]
         x = mass fraction LiBr
Outputs: h = mass specific enthalpy / [J/kg]

Based on table 7 and equation (4) in reference.
"""
    
    a=[2.27431,-7.99511, 385.239,-16394,-422.562,0.113314,-8.33474,-17383.3,\
    6.49763,3245.52,-13464.3,39932.2,-258877,-0.00193046,2.80616,-40.4479,\
    145.342,-2.74873,-449.743,-12.1794,-0.00583739,0.233910,0.341888,8.85259,\
    -17.8731,0.0735179,-0.000179430,0.00184261,-0.00624282,0.00684765]
    m=[1,1,2,3,6,1,3,5,4,5,5,6,6,1,2,2,2,5,6,7,1,1,2,2,2,3,1,1,1,1]
    n=[0,1,6,6,2,0,0,4,0,4,5,5,6,0,3,5,7,0,3,1,0,4,2,6,7,0,0,1,2,3]
    t=[0,0,0,0,0,1,1,1,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5]
    T_crit = pwater.T_critical() # [K]
    P_crit = pwater.p_critical()
    # This one had a problem starting around CoolProp version 5.0.8
    #h_crit = PropsSI('H','T',T_crit,'P',P_crit,'water') # J/kg
    # Instead, use the low-level interface
    state = AbstractState('HEOS','Water')
    state.specify_phase(constants.iphase_critical_point)
    state.update(constants.PT_INPUTS, P_crit, T_crit)
    h_c = state.hmolar() # [J/mol] = [kJ/kmol]
    T_c = T_crit # [K]
    T_0 = 221. # [K] "is a nonlinear parameter of the equations"

    TK = T
    x_N = molefraction(x)
    # print(x_N)
    s = 0
    for i in range(len(a)):
        s = s + a[i] * x_N**m[i] * (0.4-x_N)**n[i] * (T_c/(TK-T_0))**t[i]
    Qu = 0.0
    pwater.update(CP.QT_INPUTS, Qu, TK)
    h_w_molar = pwater.hmolar() - h_w_molar_ref # [J/mol]
    h_molar = (1 - x_N) * h_w_molar + h_c * s # [J/mol]
    MW = x_N * MW_LiBr + (1 - x_N) * MW_H2O # [kg/mol]
    result = h_molar / MW # [J/kg]
    return result

def massSpecificEntropy(T,x):
    """Inputs:  T = Temperature / [Kelvin]
         x = mass fraction LiBr
Outputs: s = mass specific entropy / [J/kg-K]

Based on table 8 and equation (5) in reference.
"""
    
    a=[1.53091,-4.52564, 698.302,-21666.4,-1475.33,0.0847012,-6.59523,
       -29533.1,0.00956314,-0.188679,9.31752,5.78104,13893.1,-17176.2,
       415.108,-55564.7,-0.00423409,30.5242,-1.67620,14.8283,0.00303055,
       -0.0401810,0.149252,2.59240,-0.177421,-0.0000699650,0.000605007,
       -0.00165228,0.00122966]
    m = [1,1,2,3,6,1,3,5,1,2,2,4,5,5,6,6,1,3,5,7,1,1,1,2,3,1,1,1,1]
    n = [0,1,6,6,2,0,0,4,0,0,4,0,4,5,2,5,0,4,0,1,0,2,4,7,1,0,1,2,3]
    t = [0,0,0,0,0,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5]
    T_c = 647.096 # [K]
    T_0 = 221 # [K] "is a nonlinear parameter of the equations"
    #s_c = 79.3933 # [J/gmol-K]
    T_crit = pwater.T_critical() # [K]
    P_crit = pwater.p_critical() # [Pa]
    state = AbstractState('HEOS','Water')
    state.specify_phase(constants.iphase_critical_point)
    state.update(constants.PT_INPUTS, P_crit, T_crit)
    s_crit_molar = state.smolar() # J/mol-K
    s_c = s_crit_molar
        
    TK = T
    x_N = molefraction(x)
    s = 0
    for i in range(len(a)):
         s = s + a[i] * x_N ** m[i] * (0.4 - x_N) ** n[i] \
             * (T_c / (TK - T_0)) ** t[i]
    #s_w=entropy(Water,T=T,x=0)
    Qu_water = 0.0
    #s_w_mass = PropsSI('S','T',TK,'Q',Qu_water,'water') # J/kg-K
    #s_w_molar = s_w_mass * MW_H2O # J/mol-K
    pwater.update(CP.QT_INPUTS, Qu_water, TK)
    s_w_molar = pwater.smolar() # J/mol-K
    s_molar = (1 - x_N) * s_w_molar + s_c * s
    MW = x_N * MW_LiBr + (1 - x_N) * MW_H2O
    result = s_molar / MW
    return result
    
def massSpecificHeat(T,x):
    """Inputs:  T = Temperature / [Kelvin]
         x = mass fraction LiBr
Outputs: cp = mass specific heat / [J/kg-K]

Based on Table 6 and equation (3) of reference.
"""
    a = [-14.2094,40.4943,111.135,229.980,
         1345.26,-0.0141010,0.0124977,-0.000683209]
    m = [2,3,3,3,3,2,1,1]
    n = [0,0,1,2,3,0,3,2]
    t = [0,0,0,0,0,2,3,4]
    Cp_t = 76.0226 # [J/mol-K]
    T_c=647.096 # [K]
    T_0 = 221 # [K] "is a nonlinear parameter of the equations"
    TK = T
    x_N = molefraction(x)
    s=0
    for i in range(len(a)):
        s = s + a[i] * x_N ** m[i] * (0.4 - x_N) ** n[i] \
            * (T_c / (TK - T_0)) ** t[i]
    Qu_water = 0.0
    pwater.update(CP.QT_INPUTS, Qu_water, TK)
    Cp_w_molar = pwater.cpmolar() # J/mol-K
    Cp_molar = (1 - x_N) * Cp_w_molar + Cp_t * s
    MW = x_N * MW_LiBr + (1 - x_N) * MW_H2O
    result = Cp_molar / MW
    return result
    
def twoPhaseProps(h,P,z):
    """Some notes.
    This function returns the quality, temperature and liquid composition of a
    2-phase mixture of liquid lithium bromide-water and water vapor at specific
    enthalpy h, pressure P, and overall composition, z.
    
    Inputs:
    
    h is enthalpy / [J/kg]
    P is pressure / [bar].
    z is the overall lithium bromide mass fraction [kg/kg].
    
    Outputs:
    
    T is temperature / [K].
    Q is the quality (or vapor fraction) on a mass basis [kg/kg].
    x is the lithium bromide mass fraction of the liquid phase [kg/kg].

    We observe that all the lithium bromide mass is in the liquid phase.
    Therefore, a mass balance requires (1 - Q) x = z.
    An enthalpy balance gives simply h = (1-Q) h_liquid + Q h_vapor.
    Rearranging, (h - h_liquid) = Q (h_vapor - h_liquid).
    Therefore we have two equations to solve for Q and x.
    Equilibrium requires T_liquid = T_vapor, so we can use exisiting functions.
    """
    #P_pascal = P * 1e5
    Q, T, x = 0, 0, 0
    Q = -100	# subcooled
    x = z
    T = temperature(P, x) # K
    hL = massSpecificEnthalpy(T,x) # J/kg
    if (h == hL): Q = 0
    if (h <= hL): return Q, T, x

    Q = 0.1
    for iter in range(100):
        Qlast = Q
        x = z / (1. - Q)
        T = temperature(P, x) 
        hL = massSpecificEnthalpy(T,x) # J/kg
        hv = 0
        if (h > hL):
            Q_vapor = 1.
            pwater.update(CP.QT_INPUTS, Q_vapor, T)
            hv = pwater.hmass() # J/kg
            hfg = hv - hL
            Q = (h - hL) / (hfg) # kg/kg
            # qq = (x - z) / x
        else:
            Q = 0.
        #print("{},h={},P={},z={},Q={},x={},T={},hL={},hv={}"
        #    .format(iter,h,P,z,Q,x,T,hL,hv))
        if (abs(Q - Qlast) < 0.00001) and (iter > 5):
            break
    #print "TwoPhaseProps converged at iter = ", iter
    return Q, T, x

def massSpecificGibbs(T,x):
    h = massSpecificEnthalpy(T,x) # [J/kg]
    s = massSpecificEntropy(T,x) # [J/kg-K]
    g = h - T * s # [J/kg]
    return g

def massDensity(T,x):
    """This function returns the density of a liquid lithium bromide water
solution given the temperature and composition, based on equation 2 and table
5 in Patek and Klomfar, Int. J. of Refrigeration, Vol 29, pp. 566-578, (2006).

Inputs:
    T = temperature / K
    x  = mass fraction of lithium bromide in the liquid
Outputs:
    density in units of kg/m3
"""    
    a=[1.746,4.709]
    m=[1,1]
    t=[0,6]
    
    #rho_crit_mass = PropsSI('water','rhomass_critical') # [kg/m3]
    #rho_crit_molar = rho_crit_mass / MW_H2O # [mol/m3]
    rho_crit_molar = pwater.rhomolar_critical()
    T_crit = pwater.T_critical() # [K]
    print("""By the way, water critical properties:
T_crit = {} K,
rho_crit_molar = {} mol/m3""".format(T_crit,rho_crit_molar))
    
    #rho_c = 17873 # [gmol/m^3]
    T_c=647.096 # [K]
    x_N = molefraction(x)
    # saturated liquid water density
    Qu_water = 0.0
    #rhomass_sat = PropsSI('D','T',T,'Q',Qu_water,'water') # kg/m3
    #rhomolar_sat = rhomass_sat / MW_H2O
    pwater.update(CP.QT_INPUTS, Qu_water, T)
    rhomolar_sat = pwater.rhomolar() # mol/m3

    s=0
    for i in range(len(a)):
        s = s + a[i] * (x_N ** m[i]) * ((T / T_c) ** t[i])
    d_molar = (1 - x_N) * rhomolar_sat + rho_crit_molar * s
    MW = x_N * MW_LiBr + (1 - x_N) * MW_H2O
    result = d_molar * MW # kg/m^3
    return result

# function to check whethe
crystallization_data_T = np.array(
      [ -53.6 ,  -49.32,  -42.12,  -36.32,  -32.96,  -29.17,  -25.24,
        -16.11,  -13.47,   -8.94,   -4.54,    1.11,    5.1 ,    9.93,
         18.99,   24.29,   33.14,   38.26,   44.27,   50.35,   57.58,
         63.42,   70.9 ,   71.69,   82.68,   83.11,   91.36,   91.82,
        101.05,  102.02])
crystallization_data_x = np.array(
      [ 0.452 ,  0.4803,  0.4963,  0.5009,  0.505 ,  0.512 ,  0.517 ,
        0.5195,  0.537 ,  0.5475,  0.5592,  0.5681,  0.5722,  0.5808,
        0.5867,  0.6063,  0.625 ,  0.6396,  0.6517,  0.6582,  0.6616,
        0.6655,  0.6737,  0.6739,  0.6832,  0.6827,  0.6899,  0.6905,
        0.7004,  0.7008])

x_to_T_crystal = interp1d(crystallization_data_x,crystallization_data_T)

x_min_crystal = min(x_to_T_crystal.x)
x_max_crystal = max(x_to_T_crystal.x)

def is_crystalized(T,x):
    
    if x < x_min_crystal:
        return False
    
    # a curva parece ser 'reta' por essas bandas no gráfico do tchô
    if x > x_max_crystal:
        x = x_max_crystal
    
    return T < float(x_to_T_crystal(x))


if __name__ == "__main__":
    # Unit testing
    # Pressure (cf table 6.1, or rather LiBrSS7B.EES (slight difference))
    #for TC,x in ((32.7,0.567),(63.3,0.567),(89.9,0.625),(53.3,0.625),\
    #    (44.7,62.5),(77.0,0.0),(40.2,0.0),(1.3,0.0)):
    #    P = P_LiBrH2O(C2K(TC),x)
    #    print("(T,x) = ({},{}) -> P = {}".format(TC,x,P))
    # Confer documentation for massfraction_LiBrH2O
    for x,w_expected in ((0.1718, 0.50),):
        w = mole2massFraction(x)
        print("[1] Mole fraction {} -> mass fraction {}, expected {}".format(
            x,w,w_expected))
    # Confer documentation for molefraction_LiBrH2O
    for w,x_expected in ((0.5, 0.1718),):
        x = molefraction(w)
        print("[2] Mass fraction {} -> mole fraction {}, expected {}".format(
            w,x,x_expected))
    # Confer documentation for P_LiBrH2O
    for TC,x in ((50.0,0.5),):
        T = C2K(TC)
        P = pressure(T,x)
        P_expected = 0.03486
        print("[3] (T,x) = ({} C, {}) -> P = {} bar".format(TC,x,P))
        print("Expected {}".format(P_expected))
    # Confer documentation for T_LiBrH2O
    print("[4]")
    for P,x in ((3.5e3 / 1.0e5, 0.5),):
        T = temperature(P,x)
        T_expected = 50.08
        print("[4] (P,x) = {} bar, {} -> T = {} C".format(P,x,K2C(T)))
        print("Expected {} C".format(T_expected))
        print("Checking the solution, because it is not always correct.")
        P_check = pressure(T,x)
        #print("Given T, x = {} K, {}, we get P = {} bar".format(T,x,P_check))
        #print("Do you think it failed to solve? Changing the guess...")
    
    # Confer documentation for x_LiBrH2O
    print("[5]")
    for TC,P in ((50.0, 3.5e3 / 1.0e5),):
        T = C2K(TC)
        x = massFraction(T,P)
        x_expected = 0.4995
        print("[5] (TC,P) = {}, {} -> x = {}".format(TC,P,x))
        print("Expecting {}".format(x_expected))
    # Confer documentation for h_LiBrH2O
    for TC,x,h_expected in ((50,0.5,105e3),):
        T = C2K(TC)
        h = massSpecificEnthalpy(T,x)
        print("[6] T,x = {} C, {} -> h = {} J/kg, expected {}".format(TC,x,h,h_expected))
    for TCvar in np.linspace(20,80,10):
        T = C2K(TCvar)
        h = massSpecificEnthalpy(T,x)
        print("[7] T,x = {:.3f} C, {} -> h = {} J/kg".format(TCvar,x,h))
    TC = 70.
    for xvar in np.linspace(0.2,0.8,10):
        T = C2K(TC)
        h = massSpecificEnthalpy(T,xvar)
        print("[8] T,x = {} C, {:0.3f} -> h = {} J/kg".format(TC,xvar,h))
        
    print("By the way, pure water @ 0 deg C, quality 0 has enthalpy {} J/kg"
        .format(PropsSI('H','T',273.16,'Q',0,'water')))
    T_crit = PropsSI('water','Tcrit')
    P_crit = PropsSI('water','pcrit')
    print("By the way, pure water @ critical point is has {} K, {} Pa"
        .format(T_crit, P_crit))
    state = AbstractState('HEOS','Water')
    state.specify_phase(constants.iphase_critical_point)
    state.update(constants.PT_INPUTS, P_crit, T_crit)
    h_crit = state.hmass()
    print("By the way, pure water @ critical point has enthalpy {} J/kg"
        .format(h_crit))
    
    # Confer documentation for s_LiBrH2O
    for TC,x,s_expected in ((50.0,0.5,0.3519e3),):
        T = C2K(TC)
        s = massSpecificEntropy(T,x)
        print("T,x = {} C, {} -> s = {} J/kg-K, expected {}".format(TC,x,s,s_expected))
    
    # Confer documentation for Cp_LiBrH2O
    for TC,x,Cp_expected in ((50.0,0.5,2.183e3),):
        T = C2K(TC)
        Cp = massSpecificHeat(T,x)
        print("T,x = {} C, {} -> Cp = {} J/kg-K, expected {}"
            .format(TC,x,Cp,Cp_expected))
    
    # Confer documentation for 
    P1, P2 = 10e3/1e5, 1e3/1e5 # [bar]
    TC1 = 70.0 # [C]
    x1 = 0.6 # [kg/kg]
    T1 = C2K(TC1)
    h1 = massSpecificEnthalpy(T1,x1)
    print("T={}, x={} -> h = {}, P (unused) = {}".format(T1,x1,h1,P1))
    # Throttle to lower pressure"
    h2 = h1
    z2 = x1
    QTx2 = twoPhaseProps(h2,P2,z2)
    # Solution:
    QTx2_expected = 0.0146,C2K(48.6),0.6089 # [dim], [C], [dim]
    print("Got QTx = {}, expected {}".format(QTx2, QTx2_expected))
    
    # Confer documentation for density
    TC, w = 50, 0.5
    rho, rho_expected = massDensity(C2K(TC), w), 1522 # kg/m3
    print("T, w = {} K, {} kg/kg -> rho = {} kg/m3, expected {}"
        .format(T, w, rho, rho_expected))

    ()
