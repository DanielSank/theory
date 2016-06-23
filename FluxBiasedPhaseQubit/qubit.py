from numpy import cos,sin,exp

def fp(f0,epsilon):
    return f0*(epsilon**0.5)
    
def N(f0,F0,epsilon):
    return (2.0/3)*(F0/f0)*epsilon**(5.0/2.0)
    
def f10(f0,F0,epsilon):
    return fp(f0,epsilon)*(1-((5.0/36)*(1.0/N(f0,F0,epsilon))))
    
def fnl(f0,F0,epsilon):
    return fp(f0,epsilon)*(5.0/36.0)*(1.0/N(f0,F0,epsilon))

