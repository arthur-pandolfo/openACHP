#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
# from scipy.optimize import fsolve
import matplotlib.pyplot as plt

plt.close("all")

# Patek and Klomfar, my implementation
import libr_props
myH = np.vectorize(libr_props.massSpecificEnthalpy)

# Patek and Klomfar spline fit but incorrect reference states
import CoolProp.CoolProp as CP
librname = lambda x: 'INCOMP::LiBr[{}]'.format(x)
def Hsat1(T,x):
    try:
        h = CP.PropsSI('H','T',T,'Q',0,librname(x))
        return h
    except:
        return np.nan
Hsat2 = np.vectorize(Hsat1)
def Hsat3(P,x):
    try:
        #perr = lambda(t): CP.PropsSI('P','T',t,'Q',0,librname(x))-P
        #t = fsolve(perr,360)[0]
        #h = CP.PropsSI('H','T',t,'Q',0,librname(x))
        t = libr_props.temperature(P*1e-5,x)
        h = libr_props.massSpecificEnthalpy(t,x)
        return h
    except:
        return np.nan
Hsat4 = np.vectorize(Hsat3)

# According to EES, this is from ASHRAE Fundamentals 1998.
# It is also printed in ASHRAE fundamentals 2013, chapter 30, figure 34.
# The formula is a quadratic in T and quartic in X.
# from ees_interface2 import EES_DLL
# LiBr_path = 'C:/EES32/Userlib/Libr/LIBR.dll'
# myDLL = EES_DLL(LiBr_path)
# h_libr=myDLL.func['H_LIBR']
# def hfun(T,x):
#     inarglist=[T-273.15,x*100.,2]
#     s0,outarglist=h_libr.call("",inarglist)
#     return outarglist * 1e3
# hfun2=np.vectorize(hfun)

# EES has also added recently some more stuff.
# LiBrSSC_path = 'C:/EES32/Userlib/Libr/SSCLiBr.dll'
# SSC_DLL = EES_DLL(LiBrSSC_path)
# h_librssc=SSC_DLL.func['LiBrSSCh']
# def hfun3(T,x,P=None):
#     if P==None:
#         s0,outarglist=h_librssc.call("",[T-273.15,x])
#     else:
#         s0,outarglist=h_librssc.call("",[T-273.15,x,P*1e-3])
#     return outarglist * 1e3
# hfun4=np.vectorize(hfun3)
# LiBrSSCT(P,X) ['kPa', '-'] |-> ['C']
# def t_ssc(P,x):
#     s0,outarglist = SSC_DLL.func['LiBrSSCT'].call("",[P*1e3,x])
#     return outarglist + 273.15
# t_sscv = np.vectorize(t_ssc)

x = np.linspace(0,.8)
T_ref = 293.15
h_offset = libr_props.massSpecificEnthalpy(T_ref,x)

#fig = plt.figure(figsize=(16,8))
fig = plt.figure(figsize=(8,5),dpi=600)
plt.xlabel("LiBr mass fraction, $x$ [kg/kg]")
plt.ylabel("Solution enthalpy, $h_{{solution}}$ [J/kg]")
fig.gca().ticklabel_format(axis='y', style = 'sci', scilimits=(-2,2), useOffset=False)
def _crop(h):
    if h < -1.3e5 or h > 4e5:
        return np.nan
    return h
crop = np.vectorize(_crop)
TT = [273.15,293.15,313.15,333.15,353.15,373.15]
cc = 'b g r c m y'.split()

for T,c in zip(TT,cc):
    h0 = libr_props.massSpecificEnthalpy(T,x)
    h1 = Hsat2(T,x)+h_offset
    # h2 = hfun2(T,x)
    # h4 = hfun4(T,x)
    plt.plot(x,crop(h0),color=c,label="Pátek and Klomfar" if T == 273.15 else None)
    plt.plot(x,h1,'--',label="CoolProp" if T==273.15 else None)
    # plt.plot(x,crop(h4),'--',color=c,label="Yuan and Herold" if T == 273.15 else None)
    # plt.plot(x,crop(h2),'.-',color=c,label="ASHRAE" if T == 273.15 else None)
    xtext = 0.4 + (T-273.15) * 0.001
    htext = libr_props.massSpecificEnthalpy(T+5,xtext)
    plt.annotate('{} $^\circ$C'.format(T-273.15),[xtext,htext])
bottom = 0.8 + 0*libr_props.crystallization_data_T
# h_crystal = hfun2(libr_props.crystallization_data_T+273.15,
#                 libr_props.crystallization_data_x)
# plt.fill_betweenx(h_crystal,
#                     libr_props.crystallization_data_x,
#                     bottom,
#                     edgecolor='aqua',
#                     facecolor='aqua',
#                     label="Boryta")
# plt.annotate('Crystallization (Boryta)',[0.55,-1e5])    

# plt.ylim([h_crystal.min(),4e5])
plt.xlim([0,0.8])
plt.legend(loc='upper center')

import libr3
chiller=libr3.ChillerLiBr1(T_evap=5,T_cond=45,x1=0.6026,x2=0.66)
#chiller=libr3.ChillerLiBr1(T_evap=5,T_cond=30,x1=0.5,x2=0.66)
chiller.iterate1()
chiller2=libr3.ChillerLiBr1(T_evap=5,T_cond=50,x1=0.627,x2=0.66)
chiller2.iterate1()
chiller3=libr3.ChillerLiBr1(T_evap=5,T_cond=55,x1=0.6515,x2=0.66)
chiller3.iterate1()
print(chiller)
print(chiller3)
def cyclexh(ch):
    cyclex = [ch.x1, ch.x1, ch.x2, ch.x2, ch.x1]
    cycleh = [ch.h_abs_outlet, ch.h_gen_inlet, ch.h_gen_outlet, ch.h_abs_inlet, ch.h_abs_outlet]
    return cyclex,cycleh

#fig = plt.figure(figsize=(8,5))
x=np.linspace(0.5,0.8)
fix, ax = plt.subplots(1, figsize=(8,5))
plt.xlabel("LiBr mass fraction, $x$ [kg/kg]")
plt.ylabel("Solution enthalpy, $h_{{solution}}$ [J/kg]")
ax.ticklabel_format(axis='y', style = 'sci', scilimits=(-2,2), useOffset=False)
def _crop(h):
    if h < 0 or h > 3.5e5:
        return np.nan
    return h
crop = np.vectorize(_crop)
TT = [273.15,293.15,313.15,333.15,353.15,373.15,393.15]
cc = 'b g r c m y b'.split()
for T in TT:
    h0 = libr_props.massSpecificEnthalpy(T,x)
    h1 = Hsat2(T,x)+h_offset
    ax.plot(x,crop(h0),'--',color='gray')
    xtext = 0.51
    htext = libr_props.massSpecificEnthalpy(T+5,xtext)
    ax.annotate('{} $^\circ$C'.format(T-273.15),[xtext,htext])

# Crystallization curve
bottom = 0.8 + 0*libr_props.crystallization_data_T
# h_crystal = hfun2(libr_props.crystallization_data_T+273.15,
#                 libr_props.crystallization_data_x)
# ax.fill_betweenx(h_crystal,
#                    libr_props.crystallization_data_x,
#                    bottom,
#                    edgecolor='aqua',
#                    facecolor='aqua',
#                    label="Crystallization")
#ax.annotate('Crystallization (Boryta)',[0.55,-1e5])

# Corresponds to 1 to 5 C
Trange = 2,5
Prange = [CP.PropsSI("P","T",t+273.15,"Q",0,"HEOS::water") for t in Trange]
hrange = [Hsat4(p,x) for p in Prange]
ax.fill_between(x,hrange[0],hrange[1],color='blue',alpha=0.3,
                label="Evap. pressure, {}-{}$^\circ$C".format(*Trange))
# Corresponds to 45 to 55 C
Trange = 45,55
Prange = [CP.PropsSI("P","T",t+273.15,"Q",0,"HEOS::water") for t in Trange]
hrange = [Hsat4(p,x) for p in Prange]
ax.fill_between(x,hrange[0],hrange[1],color='red',alpha=0.3,
                label="Cond. pressure, {}-{}$^\circ$C".format(*Trange))

for ch,c in zip([chiller,chiller2,chiller3],['b','g','r']):
    cyclex,cycleh = cyclexh(ch)
    ax.plot(cyclex,cycleh,'.',color=c)
    print('CHILLER',ch.x1,ch.COP)
    for i in range(4):
        ax.annotate("",(cyclex[i+1],cycleh[i+1]),(cyclex[i],cycleh[i]),
                    arrowprops=dict(arrowstyle="-|>",facecolor=c,edgecolor=c))
    ax.annotate("{}$^\circ$C".format(ch.T_cond),(ch.x1-0.013,ch.h_abs_outlet),color=c)
    # EU
    # break
ax.annotate('Absorber',(0.6,1e5))
ax.annotate('Desorber',(0.6,2.6e5))

plt.ylim([0,3.5e5])
plt.xlim([0.4,0.8])
ax.legend(loc='lower right')
# For word, export svg and convert with inkscape
#plt.savefig('../img/compare_libr_props.fig3.svg')
plt.show()