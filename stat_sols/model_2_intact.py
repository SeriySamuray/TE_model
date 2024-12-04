# expected observables under normal conditions:
a0,mq0,cq0,q0,mL0,mS0,cL0,O10,bL0,bS0,L0,S0 = np.array([1.70903634e+13, 1.01281678e+07, 3.87047959e+06, 4.75144500e+09,
                                                        1.32203563e+05, 2.01277667e+08, 2.21289804e+05, 5.75527581e+07,
                                                        4.00758209e+03, 1.16512701e+04, 1.05859770e+03, 1.32079090e+04])

# Steady state reference
sty0 = np.array([a0,mq0,cq0,q0,mL0,mS0,cL0,O10,bL0,bS0,L0,S0])

# knock-S
sty1 = np.array([4.45330154e+09, 5.01082679e+06, 2.42671386e+06, 2.97904436e+09,
                 1.03495292e+05, 0, 2.19540165e+05, 2.35291548e+05,
                 2.35291548e+05, 0, 6.18819300e+04, 0])

# knock-LS
sty2 = np.array([2.37250469e+13, 1.00338202e+07, 3.94118717e+06, 4.83824643e+09,
                 0,0,0,0,
                 0,0,0,0])

# modified: lmd_L/2
sty3 = np.array([1.09947922e+10, 7.50350728e+06, 2.89550377e+06, 3.55454517e+09,
                 6.76476378e+05, 7.41905418e+08, 1.14339623e+06, 2.99986736e+08,
                 1.01283610e+04, 6.07937141e+04, 5.34547624e+03, 6.88470715e+04])

# modified: lmd_L/2 & lmd_S/2
sty4 = np.array([6.41826251e+12, 1.01857896e+07, 3.82352127e+06, 4.69379843e+09,
                 2.21437731e+05, 7.17038794e+08, 3.64086674e+05, 9.57972649e+07,
                 2.11395419e+03, 1.93946239e+04, 1.11679887e+03, 4.39715416e+04])

# Variable parameters
A0       = 1.31e10 # energy supply [m.p.c./min]
lmd_a    = 2.87e-4 # ATP molecules degradation rate constant [1/min]
lmd_L    = 0.53    # LINE-1 deactivation rate constant [1/min]
lmd_S    = 2.47    # SINE deactivation rate constant [1/min]

# Set directories
path = "../data/"
dir_data = path+"stoch2_intact/"