import numpy
import sympy

A = sympy.Symbol('A')
B = sympy.Symbol('B')
C = sympy.Symbol('C')
D = sympy.Symbol('D')
x = sympy.Symbol('x')
m = sympy.Symbol('m')
o = sympy.Symbol('o')
a = sympy.Symbol('a')
l = sympy.Symbol('l')
k = sympy.Symbol('k')
t = sympy.Symbol('t')

class TimoshenkoBeam(object):
    def __init__(self,alpha) -> None:
        self.nu = 0.3
        self.gamma = 1/(2*(1+self.nu))*5/5
        self.alpha = alpha
        self.u_0 = None
        self.u_1 = None
        self.p_0 = None
        self.p_1 = None

    @property
    def delta(self):
        return 4*self.gamma/(1+self.gamma)**2*self.alpha/l + (1-self.gamma)**2/(1+self.gamma)**2

    @property
    def omega_squared(self):
        return 1/2*l*(1+self.gamma)*(self.delta**(1/2)+1)

    @property
    def mu_squared(self):
        return 1/2*l*(1+self.gamma)*(self.delta**(1/2)-1)

    @property
    def theta_squared(self):
        return 1/2*l*(1+self.gamma)*(1-self.delta**(1/2))

    def u(self, size_of_lam):
        if size_of_lam == 'less':
            return A*sympy.sinh(m*x) + B*sympy.cosh(m*x) + C*sympy.sinh(o*x) + D*sympy.cosh(o*x)
        elif size_of_lam == 'equal':
            return B + C*sympy.sin(o*x) + D*sympy.cos(o*x)
        elif size_of_lam == 'greater':
            return A*sympy.sin(t*x) + B*sympy.cos(t*x) + C*sympy.sin(o*x) + D*sympy.cos(o*x)

    def p(self, size_of_lam):
        if size_of_lam == 'less':
            return A*(l+m**2)/(m*sympy.cosh(m*x)) + B*(l+m**2)/(m*sympy.sinh(m*x)) + C*(-l+o**2)/(o*sympy.cosh(o*x)) + D*(l-o**2)/(o*sympy.sinh(o*x))
        elif size_of_lam == 'equal':
            return A + B*a*x + C*(-l+o**2)/(o)*sympy.cos(o*x) + D*(l-o**2)/(o)*sympy.sin(o*x)
        elif size_of_lam == 'greater':
            return A*(-l+t**2)/(t)*sympy.cos(t*x) + B*(l-t**2)/(t)*sympy.sin(t*x) + C*(-l+o**2)/(o)*sympy.cos(o*x) + D*(l-o**2)/(o)*sympy.sin(o*x)
        
    def calculate_matrix(self):
        matrix_dict = {'less': None,
                       'equal': None,
                       'greater': None}
        for key in matrix_dict:
            print(key)
            u = self.u(key)
            p = self.p(key)
            
            u_at_zero = u.subs(x, self.u_0)
            p_at_zero = p.subs(x, self.p_0)
            
            u_solved = sympy.solvers.solve(u_at_zero, [D])[0]
            p_solved = sympy.solvers.solve(p_at_zero, [C])[0]

            u = u.subs(D, u_solved)
            p = p.subs(C, p_solved)
            
            M1 = (sympy.diff(u, x, 1)).subs(x, 1)
            M2 = (sympy.diff(u, x, 1)-p).subs(x, 1)
            
            M = sympy.Matrix([[M1.subs([A, B], [1, 0]), M2.subs([A, B], [0, 1])]
                             ,[M2.subs([A, B], [1, 0]), M1.subs([A, B], [0, 1])]])
            
            L = M.subs(k, sympy.sqrt(5/6))
            L = L.subs(o, sympy.sqrt(self.omega_squared))
            L = L.subs(m, sympy.sqrt(self.mu_squared))
            L = L.subs(t, sympy.sqrt(self.theta_squared))
            
            return L
        
    def determinant(self, L):
        return L.det()
            
class CantileverTimoshenkoBeam(TimoshenkoBeam):
    def __init__(self, alpha) -> None:
        super().__init__(alpha=alpha)
        self.u_0 = 0
        self.p_0 = 0

if __name__ == "__main__":
    cantileverbeam = CantileverTimoshenkoBeam(alpha = 1200)
    L = cantileverbeam.calculate_matrix()
    function = cantileverbeam.determinant(L)
    print(function)
    