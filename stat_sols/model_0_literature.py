# expected observables under normal conditions:
a0,mq0,cq0,q0,mL0,mS0,cL0,O10,bL0,bS0,L0,S0 = np.array([5.32028879e+09, 5.43892043e+06, 2.60750605e+06, 3.20098867e+09,
                                                        5.13241818e+04, 1.21841508e+08, 1.07775213e+05, 2.78154330e+07,
                                                        2.82625278e+03, 5.64255248e+03, 1.06716935e+03, 1.33614820e+04])

# Steady state reference
sty0 = np.array([a0,mq0,cq0,q0,mL0,mS0,cL0,O10,bL0,bS0,L0,S0])

# knock-S
sty1 = np.array([2.29529555e+09, 3.36327244e+06, 1.79282544e+06, 2.20086818e+09,
                 6.99550598e+04, 0, 1.63334967e+05, 1.75569280e+05, 
                 1.75569280e+05, 0, 6.59470834e+04, 0])

# knock-LS
sty2 = np.array([5.96205458e+09,  5.69819032e+06,  2.72706015e+06,  3.34775568e+09,
                 0,0,0,0,
                 0,0,0,0])

# modified: lmd_L/2
sty3 = np.array([3.60500385e+09, 4.52734948e+06, 2.18479470e+06, 2.68205969e+09,
                 2.29815368e+05, 4.55392447e+08, 4.85769072e+05, 1.27407054e+08,
                 4.47819184e+03, 2.58724875e+04, 3.37858095e+03, 6.12055856e+04])

# modified: lmd_L/2 & lmd_S/2
sty4 = np.array([4.27240077e+09, 4.89876923e+06, 2.39235807e+06, 2.93686833e+09,
                 9.18023587e+04, 4.30559873e+08, 1.96370953e+05, 5.15602878e+07,
                 1.58108214e+03, 1.04653946e+04, 1.19341160e+03, 4.95387336e+04])

# Variable parameters
A0       = 1.31e10 # energy supply [m.p.c./min]
lmd_a    = 1.47    # ATP molecules degradation rate constant [1/min]
lmd_L    = 0.37    # LINE-1 deactivation rate constant [1/min]
lmd_S    = 1.18    # SINE deactivation rate constant [1/min]

# Set directories
path = "../data/"
dir_data = path+"stoch/"