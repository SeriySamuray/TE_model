# expected observables under normal conditions:
a0,mq0,cq0,q0,mL0,mS0,cL0,O10,bL0,bS0,L0,S0 = np.array([1.98770116e+13, 1.01285391e+07, 3.87055894e+06, 4.75154242e+09,
                                                        1.32210482e+05, 2.01292320e+08, 2.21297810e+05, 5.75549012e+07,
                                                        4.00747957e+03, 1.16517030e+04, 1.05857071e+03, 1.32084009e+04])

# Steady state reference
sty0 = np.array([a0,mq0,cq0,q0,mL0,mS0,cL0,O10,bL0,bS0,L0,S0])

# knock-S
sty1 = np.array([4.45337511e+09, 5.01542216e+06, 2.42335112e+06, 2.97491625e+09,
                 1.12773512e+05, 0, 2.38671282e+05, 2.55756653e+05,
                 2.55756653e+05, 0, 6.72744717e+04, 0])

# knock-LS
sty2 = np.array([2.65122462e+13,  1.00340175e+07,  3.94123250e+06,  4.83830209e+09,
                 0,0,0,0,
                 0,0,0,0])

# modified: lmd_L/2
sty3 = np.array([1.33664548e+10, 7.93398062e+06, 2.98020632e+06, 3.65852828e+09,
                 7.31734616e+05, 8.26330148e+08, 1.20390745e+06, 3.16024585e+08,
                 1.00070113e+04, 6.40326382e+04, 5.28236831e+03, 7.25279187e+04])

# modified: lmd_L/2 & lmd_S/2
sty4 = np.array([9.19942271e+12, 1.01879215e+07, 3.82397573e+06, 4.69435633e+09,
                 2.21473515e+05, 7.17227417e+08, 3.64112583e+05, 9.58041946e+07,
                 2.11364759e+03, 1.93960169e+04, 1.11663747e+03, 4.39747226e+04])

# Variable parameters
A0       = 1.39e10 # energy supply [m.p.c./min]
lmd_a    = 2.87e-4 # ATP molecules degradation rate constant [1/min]
lmd_L    = 0.53    # LINE-1 deactivation rate constant [1/min]
lmd_S    = 2.47    # SINE deactivation rate constant [1/min]

# Set directories
path = "../data/"
dir_data = path+"stoch2_intact_azaC/"