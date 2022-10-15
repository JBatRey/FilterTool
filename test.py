from cmath import pi
from backend import *
import numpy as np

""" 

    Tipos de filtro:
    
    'butter'
    'cheby'
    'cheby2'
    'cauer'
"""

"""
filter = 'cheby2'
filter_type = 'lowpass'
Wpass = 150
Watt = 300
Gp = -3
Ga = -50
deg = 0.5
"""
"""
filter = 'cheby2'
filter_type = 'highpass'
Watt = 150
Wpass = 300
Gp = -3
Ga = -50
deg = 0.5
"""
"""
filter = 'cheby2'
filter_type = 'bandpass'
Watt = [100,350]
Wpass = [10,300]
Gp = -3
Ga = -50
deg = 0.5
"""

filter = "cheby2"
filter_type = "bandstop"
Wpass = [2 * np.pi * (37e3 - 5e3), 2 * np.pi * (37e3 + 5e3)]
Watt = [2 * np.pi * (37e3 - 300), 2 * np.pi * (37e3 + 300)]
Gp = -6
Ga = -50
den = 0

order, Wn = get_min_order(filter, Wpass, Watt, -Gp, -Ga)

b, a = get_filter(filter, filter_type, order + 1, Wn, Wpass, Watt, Gp, Ga, den)
return_p_z(b, a)

"""Hay que elegir entre el grafico d etoda la vida, o de atenuación (al revés)
los parametros son 'attenuation' y 'standard' """
graph_filter("standard", filter, filter_type, order, Wpass, Watt, Gp, Ga, 0.5, b, a)
print(signal.TransferFunction(b, a))
