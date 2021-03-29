#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import libr3

from scipy.optimize import fsolve

### exemplo

# x1 e x2 concentracoes do soluto - dao as temperaturas do absorvedor e condensador
# m_pump é a vazão da bomba de soluto, que altera os calores trocados no evaporador
# e condensador - deixando 1 o resultado em W e W/kg é o mesmo

chiller=libr3.ChillerLiBr1(T_evap=5,T_cond=45,x1=0.6026,x2=0.66,m_pump=1.0,Eff_SHX=0.64)
chiller.iterate1()

# aqui tem todos os resultados numa tabela, e cada entrada da tabela pode
# ser conseguida como 'chiller.entrada'
# print(chiller)

# não vi ainda se tem como se fixar essas temperaturas, mas daria para iterar
# as concentrações para se conseguir
print("T gerador (inlet) [C]",chiller.T_gen_inlet)
print("T absorvedor (outlet) [C]",chiller.T_abs_outlet_max)
print("COP", chiller.COP)

# digamos que o calor a ser removido do evaporador seja 0.5 kg/s leite
Q_leite = 0.5*3900*(40-5)
print('m_dot bomba soluto', Q_leite/chiller.Q_evap_heat)

### exemplo 2: fixando-se as temperaturas do condensador e evaporador
# nao estava implementado no código original (confere chiller.iterate2())
# a ideia é obter as concentrações dados P e T para usá-las no código existente

import libr_props
from CoolProp.CoolProp import PropsSI
from scipy.interpolate import interp1d

# o que acho que seria o método implementado para achar x dado P,T
# (libr_props.massFraction) parece que não retorna os valores originais, então
# esta função get_x faria o serviço

# se o fsolve der problema de convergência o resultado está errado (não dá erro,
# só printa um warning)

def get_x(T,P):
    # T em C, P em Pa
    P_bar = P*1e-5
    T_K = libr_props.C2K(T)
    def x_res(x):
        return T_K-libr_props.temperature(P_bar,x)
    x = float(fsolve(x_res,0.6))
    return x

# crystalização
x_to_T_crystal = interp1d(libr_props.crystallization_data_x,
                          libr_props.crystallization_data_T)

x_min_crystal = min(x_to_T_crystal.x)
x_max_crystal = max(x_to_T_crystal.x)

def is_crystalized(T,x):
    
    if x < x_min_crystal:
        return False
    
    # a curva parece ser 'reta' por essas bandas no gráfico do tchô
    if x > x_max_crystal:
        x = x_max_crystal
    
    return T < float(x_to_T_crystal(x))

def check_crystal(chiller):
    x1,x2 = chiller.x1,chiller.x2
    pts = (
        ('T_SHX_concentrate_outlet',x2),
        ('T_abs_inlet_max',x2),
        ('T_abs_outlet_max',x1),
        #('T_abs_pre',),
        #('T_cond',),
        #('T_evap',),
        ('T_gen_inlet',x1),
        ('T_gen_outlet',x2))
        #('T_gen_pre',))
        
    is_ok = True
        
    for T_str,x in pts:
        
        T = getattr(chiller,T_str)
        
        if is_crystalized(T,x):
            print(T_str,('%f,%s'%(T,x)),'lies at crystallization')
            is_ok = False
            
    return is_ok

# valores
T_abs_outlet_max,T_gen_outlet = 45.,120.
T_evap,T_cond = 5,45

P_evap,P_cond = PropsSI('P',
                        'T',[libr_props.C2K(T_evap),libr_props.C2K(T_cond)],
                        'Q',[1.,1.],
                        'Water')

x1,x2 = get_x(T_abs_outlet_max,P_evap),get_x(T_gen_outlet,P_cond)

chiller=libr3.ChillerLiBr1(
    T_evap=T_evap,
    T_cond=T_cond,
    x1=x1,
    x2=x2,
    m_pump=1.0,
    Eff_SHX=0.64)
chiller.iterate1()

# ZeroCheck confere se a energia foi conservada no sistema
if abs(chiller.ZeroCheck()) > 1.0:
    raise ValueError("Parece que divergiu")
if not check_crystal(chiller):
    raise ValueError("Crystallization has happened!")