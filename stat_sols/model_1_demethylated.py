# expected observables under normal conditions:
a0,mq0,cq0,q0,mL0,mS0,cL0,O10,bL0,bS0,L0,S0 = np.array([2.40931696e+12, 1.01120553e+07, 3.86703330e+06, 4.74721425e+09,
                                                        1.31903161e+05, 2.00641730e+08, 2.20941873e+05, 5.74596157e+07,
                                                        4.01203844e+03, 1.16324573e+04, 1.05977067e+03, 1.31865310e+04])

# Steady state reference
sty0 = np.array([a0,mq0,cq0,q0,mL0,mS0,cL0,O10,bL0,bS0,L0,S0])

# knock-S
sty1 = np.array([4.45330083e+09, 5.01078199e+06, 2.42674665e+06, 2.97908462e+09,
                 1.03404963e+05, 0, 2.19353478e+05, 2.35091844e+05, 
                 2.35091844e+05, 0, 6.18293081e+04, 0])

# knock-LS
sty2 = np.array([5.77673751e+09, 9.84214464e+04, 6.79595113e+06, 8.34274993e+09,
                 0,0,0,0,
                 0,0,0,0])

# modified: lmd_L/2
sty3 = np.array([1.09450724e+10, 7.49307879e+06, 2.89337860e+06, 3.55193626e+09,
                 6.75168280e+05, 7.39911436e+08, 1.14193475e+06, 2.99599256e+08,
                 1.01318289e+04, 6.07154616e+04, 5.34728221e+03, 6.87581407e+04])

# modified: lmd_L/2 & lmd_S/2
sty4 = np.array([9.20408488e+11, 1.01438397e+07, 3.81456247e+06, 4.68280039e+09,
                 2.20733025e+05, 7.13326458e+08, 3.63575006e+05, 9.56604101e+07,
                 2.12001272e+03, 1.93671134e+04, 1.11998811e+03, 4.39087201e+04])

# Variable parameters
A0       = 1.31e10 # energy supply [m.p.c./min]
lmd_a    = 2.04e-3 # ATP molecules degradation rate constant [1/min]
lmd_L    = 0.53    # LINE-1 deactivation rate constant [1/min]
lmd_S    = 2.47    # SINE deactivation rate constant [1/min]

# Set directories
path = "../data/"
dir_data = path+"stoch1/"