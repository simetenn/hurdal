from __future__ import division
from numpy import *
from pylab import *
import chaospy as cp

## Functions
# K channel
alpha_n = vectorize(lambda v: 0.01*(-v + 10)/(exp((-v + 10)/10) - 1) if v != 10 else 0.1)
beta_n  = lambda v: 0.125*exp(-v/80)
n_inf   = lambda v: alpha_n(v)/(alpha_n(v) + beta_n(v))

# Na channel (activating)
alpha_m = vectorize(lambda v: 0.1*(-v + 25)/(exp((-v + 25)/10) - 1) if v != 25 else 1)
beta_m  = lambda v: 4*exp(-v/18)
m_inf   = lambda v: alpha_m(v)/(alpha_m(v) + beta_m(v))

# Na channel (inactivating)
alpha_h = lambda v: 0.07*exp(-v/20)
beta_h  = lambda v: 1/(exp((-v + 30)/10) + 1)
h_inf   = lambda v: alpha_h(v)/(alpha_h(v) + beta_h(v))

### channel activity ###
v = arange(-50,151) # mV


## setup parameters and state variables
T     = 45    # ms
dt    = 0.025 # ms
time  = arange(0,T+dt,dt)

## HH Parameters
V_rest  = 0      # mV
Cm      = 1      # uF/cm2
gbar_Na = 120    # mS/cm2
gbar_K  = 36     # mS/cm2
gbar_l  = 0.3    # mS/cm2
E_Na    = 115    # mV
E_K     = -12    # mV
E_l     = 10.613 # mV

Vm      = zeros(len(time)) # mV
Vm[0]   = V_rest
m       = m_inf(V_rest)      
h       = h_inf(V_rest)
n       = n_inf(V_rest)

## Stimulus
I = zeros(len(time))
for i, t in enumerate(time):
  if 5 <= t <= 30: I[i] = 10 # uA/cm2

## Simulate Model

def u(gbar_Na, gbar_K, gbar_l):
    Vm      = zeros(len(time)) # mV
    Vm[0]   = V_rest
    m       = m_inf(V_rest)      
    h       = h_inf(V_rest)
    n       = n_inf(V_rest)
    
    for i in range(1,len(time)):
        g_Na = gbar_Na*(m**3)*h
        g_K  = gbar_K*(n**4)
        g_l  = gbar_l

        m += dt*(alpha_m(Vm[i-1])*(1 - m) - beta_m(Vm[i-1])*m)
        h += dt*(alpha_h(Vm[i-1])*(1 - h) - beta_h(Vm[i-1])*h)
        n += dt*(alpha_n(Vm[i-1])*(1 - n) - beta_n(Vm[i-1])*n)

        Vm[i] = Vm[i-1] + (I[i-1] - g_Na*(Vm[i-1] - E_Na) - g_K*(Vm[i-1] - E_K) - g_l*(Vm[i-1] - E_l)) / Cm * dt
    return Vm


gbar_Na = cp.Normal(120,12)
gbar_K  = cp.Normal(36,03.6)
gbar_l  = cp.Normal(0.3,0.003)
dist = cp.J(gbar_Na, gbar_K, gbar_l)

P = cp.orth_ttr(3, dist)
nodes = dist.sample(2*len(P), "M")
solves = [u(*s) for s in nodes.T]
U_hat = cp.fit_regression(P, nodes, solves,rule="T")


s = dist.sample(10**3)
u_mc = U_hat(*s)
mean = np.mean(u_mc,1)
var = np.var(u_mc,1)



    
p_10 = np.percentile(u_mc,10,1)
p_90 = np.percentile(u_mc,90,1)

#figure()
#plot(time, var)
#plot(time, cp.Var(U_hat,dist))

figure()
plot(time, u(120*0.5,36*0.5,0.3*0.5), "y",linewidth=2 )
plot(time, u(120*1.5,36*1.5,0.3*1.5), "b",linewidth=2 )
plot(time, u(120,36,0.3), "r",linewidth=2 )
savefig("hh2.pdf")
figure()
plot(time, u(120,36,0.3), "r",linewidth=2 )
savefig("hh1.pdf")

figure()

## plot membrane potential trace
plot(time, u(120,36,0.3), "y",linewidth=2 )
plot(time, cp.E(U_hat,dist),"b",linewidth=2)
#plot(time, cp.Var(U_hat,dist))

plot(time, p_10,"r",linewidth=1)
fill_between(time, p_10,p_90,alpha=0.25)
plot(time, p_90,"g",linewidth=1)


title('Hodgkin-Huxley model for the action potential')
ylabel('Membrane Potential [mV]')
xlabel('Time [ms]')
xlim([0,T])
#ylim([-35,110])
#legend(["Known parameters", "Uncertain parameters","10 percentile", "90 percentile"])
#rc("figure",figsize=[6,4])
savefig("potential.png")


S_Ti = cp.Sens_t(U_hat, dist)

figure()
plot(time, S_Ti[0],linewidth=2)
plot(time, S_Ti[1],linewidth=2)
plot(time, S_Ti[2],linewidth=2)
xlabel("Time")
ylabel("Sensitivity")
legend(["gbar_Na","gbar_K","gbar_l"])

show()
