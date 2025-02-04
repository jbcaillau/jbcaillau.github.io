import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.linalg import lu_factor, lu_solve

def hamilton_system(t, x):
    xx, xy, xz , px, py, pz = x
    H = np.sqrt((pz * xy - py * xz)**2 + (px * xz - pz * xx)**2)
    if H == 0 : return x
    ux = (pz * xy - py * xz) / H
    uy = (px * xz - pz * xx) / H

    dx = uy * xz
    dy = -ux * xz
    dz = -uy * xx + ux * xy
    dpx = uy * pz
    dpy = -ux * pz
    dpz = ux *py - uy * px 
    return np.array([dx, dy, dz, dpx, dpy, dpz])



def shooting_system(p0_t):
    s = solve_ivp(hamilton_system, [t0, p0_t[-1]],np.concatenate((x0, p0_t[:3])) , method='RK45', t_eval=np.linspace(t0, p0_t[-1], 500))
    leq = np.linalg.norm(p0_t[:3])**2 -1
    return [s.y[0,-1]-xf[0],s.y[1,-1]-xf[1],s.y[2,-1]-xf[2],leq]

# Font: https://math.stackexchange.com/questions/2869047/python-compute-jacobian-numerically
def jacobian(F, x, delta=1e-6):
    n = len(x)
    J = np.zeros((n, n))
    for i in range(n):
        xp = np.array(x, dtype=float)
        xm = np.array(x, dtype=float)
        xp[i] += delta
        xm[i] -= delta
        J[:, i] = (np.array(F(xp)) - np.array(F(xm))) / (2 * delta)
    return J

def newton_method(F,p0t_guess, tol=1e-3, max_iter=100):
    p0_t = np.array(p0t_guess,dtype = float)  
    
    for iteration in range(max_iter):
        S_val = np.array(F(p0_t))  
        if np.linalg.norm(S_val) < tol:
            return p0_t
        
        J = jacobian(F, p0_t)

        try:
            lu, piv = lu_factor(J)  
            delta = lu_solve((lu, piv), -S_val)
        except np.linalg.LinAlgError as e:
            raise RuntimeError(f"Error: {e}")
        
        
        # p0 = np.array([delta[2],delta[1],delta[0]])
        # p0 /= np.linalg.norm(p0)
        # t = delta[3]
        # p0_t = np.array(p0.tolist() + [t])
        p0_t += delta
        # Normalizar p0
        p0_t[:3] /= np.linalg.norm(p0_t[:3])
        # print(f"I {iteration}, ||S|| = {np.linalg.norm(S_val):.2f}, t: {p0_t[-1]}, p0: {p0_t[:3]}")
    return p0_t



#information of the problem
x0 = [1.0, 0.0, 0.0]  
xf = [0.0, 1.0, 0.0]
p0 = [0.0,1/np.sqrt(3),1.0]
t0 = 0.0  
tf = np.pi * np.sqrt(3)/2 
h = 0.01
initial =  list(np.concatenate((x0, p0))) 

#initial configuration
p0_initial = [0.0,1/np.sqrt(3),1.0]
tf_initial = np.pi * np.sqrt(3)/2 

p0_initial = [1000., -0.33 , -0.57]
p0_initial = [0.74, 0.33 , 0.57]
tf_initial = 3.


#Apply newton method 
# p0_initial /= np.linalg.norm(p0_initial)
p0_t_initial = p0_initial + [tf_initial]
sol = newton_method(shooting_system, p0_t_initial)

print("Optimal solution:")
print("p0 =", sol[:3])
print("tf =", sol[-1])



solution_corrected = solve_ivp(hamilton_system, [t0, sol[-1]], np.concatenate((x0, sol[:3])), method='RK45', t_eval=np.linspace(t0, sol[-1], 500))

t_corrected = solution_corrected.t
x_corrected = solution_corrected.y.T

plt.figure("Shooting Algorithm")
for i in range(3):
    plt.plot(t_corrected, x_corrected[:, i], label=f'x{i+1}(t)')
plt.xlabel('Time (t)')
plt.ylabel('x(t)')
plt.legend()
plt.title("Shooting Algorithm")
plt.grid()
plt.show()



solution = solve_ivp(hamilton_system, [t0, tf], initial, method='RK45', t_eval=np.linspace(t0, tf, 500))

t_scipy = solution.t
x_scipy = solution.y.T

plt.figure("RK45")
for i in range(3):
    plt.plot(t_scipy, x_scipy[:, i], label=f'x{i+1}(t)')
plt.xlabel('Time (t)')
plt.ylabel('x(t)')
plt.title("RK45 with optimal initial solutions")
plt.legend()
plt.grid()
plt.show()


