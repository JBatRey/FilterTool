from backend import *



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

filter = 'cauer'
filter_type = 'bandstop'
Wpass = [100,350]
Watt = [150,300]
Gp = -3
Ga = -50
deg = 0.5


order, Wn = get_min_order(filter, Wpass, Watt, -Gp, -Ga)

b, a = get_filter(filter,filter_type, order, Wn,  Wpass, Watt, Gp,Ga,0.5)
return_p_z(b,a)
graph_filter('attenuation',filter,filter_type, order, Wpass, Watt, Gp,Ga,0.5,b,a)
print(signal.TransferFunction(b,a))

