# TE_model
This repository contains code with a model describing the basic biochemical processes in a cell in conditions of limited bioenergetic resources and competition for ATP and ribosomes from retrotransposons.
The model separately considers the processes associated with housekeeping genes and two types of non-LTR retrotransposons: LINEs and SINEs.

## Language and description of contents
We use Jupiter notebooks for interactive charting and calculations.
Up to date **(04.12.2024)** this repository contains following 5 notebooks:
1. **TE_model_det.ipynb** - for a preliminary study of stationary solutions (finding roots).
2. **TE_model_symbolic.ipynb** - study of the stability of the found stationary solutions.
3. **TE_model_stochastic.ipynb** - study the stochastic version of the model and plot graphs.
4. **TE_model_calc.ipynb** - main notebook for plotting graphs.
5. **TE_model_compare.ipynb** - compare cell states for reference and modified (assuming demethylation) models.

Folder **stat_sols/** contains data with changable model parameters and stationary solutions.
These are **.py** scripts (one script per model), executed from Jupiter notebooks.

## Description of the mathematical model

### Dynamical variables:
```math
\begin{flalign*}
\boldsymbol{a}  &\text{- amount of available energy (ATP)} && \\
\boldsymbol{mq} &\text{- the RNA of housekeeping genes} && \\
\boldsymbol{cq} &\text{- complexes of \textbf{mq} with ribosomes} && \\
\boldsymbol{q}  &\text{- proteins translated from \textbf{cq}} && \\
\boldsymbol{mL} &\text{- RNA of LINE-1} && \\
\boldsymbol{mS} &\text{- RNA of SINE} && \\
\boldsymbol{cL} &\text{- complexes of \textbf{mL} with ribosomes} && \\
\boldsymbol{O_1}&\text{- ORF1p proteins translated from \textbf{cL} (needed for new LINE-1 integration)} && \\
\boldsymbol{bL} &\text{- complexes of ORF2p with \textbf{mL}} && \\
\boldsymbol{bS} &\text{- complexes of ORF2p proteins with \textbf{mS}} && \\
\boldsymbol{L}  &\text{- LINE-1 transposons} && \\
\boldsymbol{S}  &\text{- SINE transposons} && \\
\end{flalign*}
```

### Ordinary Differential Equations
[comment]: <> (equations 1-4)
```math
\begin{flalign*}
& 1. \quad \frac{\boldsymbol{da}}{\boldsymbol{dt}} = A_0 - \lambda_{a}\boldsymbol{a} -
   v_{repl}(\boldsymbol{a},\boldsymbol{bL},\boldsymbol{O_1},\boldsymbol{bS}) -
   N_{nt} \left(N_{Q}\omega_{q}(\boldsymbol{a}) + N_{L}\omega_{L}(\boldsymbol{L},\boldsymbol{a}) + N_{S}\omega_{S}(\boldsymbol{S},\boldsymbol{a})\right) - && \\
   & \qquad \qquad N_{aa} \left(N_{q}v_{q}(\boldsymbol{cq},\boldsymbol{a}) + \frac{N_L}{3}v_L(\boldsymbol{cL},\boldsymbol{a})\right) -
   N_{nt} \left(N_{L}v_{int_L}(\boldsymbol{bL},\boldsymbol{O_1},\boldsymbol{a}) + N_{S}v_{int_S}(\boldsymbol{bL},\boldsymbol{a})\right) && \\

& 2. \quad \frac{\boldsymbol{dmq}}{\boldsymbol{dt}} = \omega_{q}(\boldsymbol{a}) - k_{bq}\:f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\:\boldsymbol{mq} +
   v_{q}(\boldsymbol{cq},\boldsymbol{a}) + k_{uq}\boldsymbol{cq} - d_{mq}\boldsymbol{mq} && \\

& 3. \quad \frac{\boldsymbol{dcq}}{\boldsymbol{dt}} = k_{bq}\:f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\:mq - k_{uq}\boldsymbol{cq} - v_{q}(\boldsymbol{cq},\boldsymbol{a}) - d_{cq}\boldsymbol{cq} && \\

& 4. \quad \frac{\boldsymbol{dq}}{\boldsymbol{dt}} = v_{q}(\boldsymbol{cq},\boldsymbol{a}) - d_{q}\boldsymbol{q} &&
\end{flalign*}
```

[comment]: <> (equations 5-8)
```math
\begin{flalign*}
& 5. \quad \frac{\boldsymbol{dmL}}{\boldsymbol{dt}} = \omega_{L}(\boldsymbol{L},\boldsymbol{a}) - k_{bL}f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\:\boldsymbol{mL} +
   k_{uL}\boldsymbol{cL} + k_{sub_S}\boldsymbol{mS}\:\boldsymbol{bL} - k_{sub_L}\boldsymbol{mL}\:\boldsymbol{bS} - d_{mL}\boldsymbol{mL} && \\

& 6. \quad \frac{\boldsymbol{dmS}}{\boldsymbol{dt}} = \omega_{S}(\boldsymbol{S},\boldsymbol{a}) - k_{sub_S}\boldsymbol{mS}\:\boldsymbol{bL} + k_{sub_L}\boldsymbol{mL}\:\boldsymbol{bS} - d_{mS}\boldsymbol{mS} && \\

& 7. \quad \frac{\boldsymbol{dcL}}{\boldsymbol{dt}} = k_{bL}f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\boldsymbol{mL} - k_{uL}\boldsymbol{cL} - v_{L}(\boldsymbol{cL},\boldsymbol{a}) - d_{cL}\boldsymbol{cL} && \\

& 8. \quad \frac{\boldsymbol{dO1}}{\boldsymbol{dt}} = v_{L}(cL,a) - v_{int_L}(\boldsymbol{bL},\boldsymbol{01},\boldsymbol{a}) - d_{O1}\boldsymbol{O1} &&
\end{flalign*}
```

[comment]: <> (equations 9-12)
```math
\begin{flalign*}
& 9. \quad \frac{\boldsymbol{dbL}}{\boldsymbol{dt}} = v_{L}(cL,a) - v_{int_L}(\boldsymbol{bL},\boldsymbol{O1},\boldsymbol{a}) - k_{sub_S}\boldsymbol{mS}\:\boldsymbol{bL} +
    k_{sub_L}\boldsymbol{mL}\:\boldsymbol{bS} - d_{bL}\boldsymbol{bL} && \\

& 10.\quad \frac{\boldsymbol{dbS}}{\boldsymbol{dt}} = k_{sub_S}\boldsymbol{mS}\:\boldsymbol{bL} - k_{sub_L}\boldsymbol{mL}\:\boldsymbol{bS} -
    v_{int_S}(\boldsymbol{bS},\boldsymbol{a}) - d_{bS}\boldsymbol{bS} && \\

& 11.\quad \frac{\boldsymbol{dL}}{\boldsymbol{dt}} = v_{int_L}(\boldsymbol{bL},\boldsymbol{01},\boldsymbol{a}) - \lambda_{L}\boldsymbol{L} && \\

& 12.\quad \frac{\boldsymbol{dS}}{\boldsymbol{dt}} = v_{int_S}(\boldsymbol{bS},\boldsymbol{a}) - \lambda_{S}\boldsymbol{S} && \\
\end{flalign*}
```

### Processes with propensities
| Process                                           | Reaction                                                                                                                      | Propensity                                                                    |
| :---                                              | :---                                                                                                                          | :---                                                                          |
| Energy supply                                     | $\emptyset\xrightarrow{A_0}\boldsymbol{a}\quad$                                                                               | $A_0$                                                                         |
| Replication                                       | $\boldsymbol{a}\xrightarrow{v_{repl}(\boldsymbol{a},\boldsymbol{bL},\boldsymbol{O_1},\boldsymbol{bS})}\emptyset\quad$         | $v_{repl}(\boldsymbol{a},\boldsymbol{bL},\boldsymbol{O_1},\boldsymbol{bS})$   |
| Transcription of housekeeping genes               | $\emptyset\xrightarrow{\omega_q(\boldsymbol{a})}\boldsymbol{mq}\quad$                                                         | $\omega_{q}(\boldsymbol{a})$                                                  |
| Transcription of LINE-1                           | $\emptyset\xrightarrow{\omega_L(\boldsymbol{L},\boldsymbol{a})}\boldsymbol{mL}\quad$                                          | $\omega_{L}(\boldsymbol{L},\boldsymbol{a})$                                   |
| Transcription of SINE                             | $\emptyset\xrightarrow{\omega_S(\boldsymbol{S},\boldsymbol{a})}\boldsymbol{mS}\quad$                                          | $\omega_{S}(\boldsymbol{S},\boldsymbol{a})$                                   |
| Ribosome binding with housekeeping genes RNA      | $\boldsymbol{mq}+f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\xrightarrow{k_{bq}}\boldsymbol{cq}\quad$                            | $k_{bq}f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\boldsymbol{mq}$               |
| Ribosome binding with LINE-1 RNA                  | $\boldsymbol{mL}+f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\xrightarrow{k_{bL}}\boldsymbol{cL}\quad$                            | $k_{bL}f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\boldsymbol{mL}$               |
| Ribosome unbinding housekeeping genes RNA         | $\boldsymbol{mq}+f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\xleftarrow{k_{uq}}\boldsymbol{cq}\quad$                             | $k_{uq}\boldsymbol{cq}$                                                       |
| Ribosome unbinding LINE-1 RNA                     | $\boldsymbol{mL}+f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\xleftarrow{k_{uL}}\boldsymbol{cL}\quad$                             | $k_{uL}\boldsymbol{cL}$                                                       |
| Translation of housekeeping genes RNAs            | $\boldsymbol{cq}\xrightarrow{v_q(\boldsymbol{cq},\boldsymbol{a})}\boldsymbol{q}\quad$                                         | $v_{q}(\boldsymbol{cq},\boldsymbol{a})$                                       |
| Translation of LINE-1 RNAs                        | $\boldsymbol{cL}\xrightarrow{v_L(\boldsymbol{cL},\boldsymbol{a})}\boldsymbol{bL}+\boldsymbol{O_1}\quad$                       | $v_{L}(\boldsymbol{cL},\boldsymbol{a})$                                       |
| Replacement $mL$ --> $mS$ in complex with $O_2$   | $\boldsymbol{mS}+\boldsymbol{bL}\xrightarrow{k_{sub_S}}\boldsymbol{bS}+\boldsymbol{mL}\quad$                                  | $k_{sub_S}\boldsymbol{mS}\boldsymbol{bL}$                                     |
| Replacement $mS$ --> $mL$ in complex with $O_2$   | $\boldsymbol{mL}+\boldsymbol{bS}\xrightarrow{k_{sub_L}}\boldsymbol{bL}+\boldsymbol{mS}\quad$                                  | $k_{sub_L}\boldsymbol{mL}\boldsymbol{bS}$                                     |
| Integration of LINE-1                             | $\boldsymbol{bL}+\boldsymbol{0_1}\xrightarrow{v_{int_L}(\boldsymbol{bL},\boldsymbol{01},\boldsymbol{a})}\boldsymbol{L}\quad$  | $v_{int_L}(\boldsymbol{bL},\boldsymbol{01},\boldsymbol{a})$                   |
| Integration of SINE                               | $\boldsymbol{bS}\xrightarrow{v_{int_S}(\boldsymbol{bS},\boldsymbol{a})}\boldsymbol{S}\quad$                                   | $v_{int_S}(\boldsymbol{bS},\boldsymbol{a})$                                   |
| ATP degradation                                   | $\boldsymbol{a} \xrightarrow{\lambda_{a}}\emptyset\quad$                                                                      | $\lambda_{a}\boldsymbol{a}$                                                   |
| $mq$ degradation                                  | $\boldsymbol{bS}\xrightarrow{d_{mq}}\emptyset\quad$                                                                           | $d_{mq}\boldsymbol{mq}$                                                       |
| $cq$ degradation                                  | $\boldsymbol{bL}\xrightarrow{d_{cq}}\emptyset\quad$                                                                           | $d_{cq}\boldsymbol{cq}$                                                       |
| housekeeping genes degradation                    | $\boldsymbol{O1}\xrightarrow{d_{q} }\emptyset\quad$                                                                           | $d_{q}\boldsymbol{q}$                                                         |
| $mL$ degradation                                  | $\boldsymbol{cL}\xrightarrow{d_{mL}}\emptyset\quad$                                                                           | $d_{mL}\boldsymbol{mL}$                                                       |
| $mS$ degradation                                  | $\boldsymbol{mS}\xrightarrow{d_{mS}}\emptyset\quad$                                                                           | $d_{mS}\boldsymbol{mS}$                                                       |
| $cL$-complex degradation                          | $\boldsymbol{mL}\xrightarrow{d_{cL}}\emptyset\quad$                                                                           | $d_{cL}\boldsymbol{cL}$                                                       |
| $O_1$ degradation                                 | $\boldsymbol{q} \xrightarrow{d_{O1}}\emptyset\quad$                                                                           | $d_{O1}\boldsymbol{O1}$                                                       |
| $bL$-complex degradation                          | $\boldsymbol{cq}\xrightarrow{d_{bL}}\emptyset\quad$                                                                           | $d_{bL}\boldsymbol{bL}$                                                       |
| $bS$-complex degradation                          | $\boldsymbol{mq}\xrightarrow{d_{bS}}\emptyset\quad$                                                                           | $d_{bS}\boldsymbol{bS}$                                                       |
| LINE-1 degradation                                | $\boldsymbol{L} \xrightarrow{\lambda_{L}}\emptyset\quad$                                                                      | $\lambda_{L}\boldsymbol{L}$                                                   |
| SINE degradation                                  | $\boldsymbol{S} \xrightarrow{\lambda_{S}}\emptyset\quad$                                                                      | $\lambda_{S}\boldsymbol{S}$                                                   |

### Some reaction rates

#### Replication
$\quad v_{repl}(\boldsymbol{a},\boldsymbol{bL},\boldsymbol{O_1},\boldsymbol{bS}) = N_{nt}\left(\frac{N_{g}}{\tau} + N_{L}v_{int_L}(\boldsymbol{bL},\boldsymbol{O_1},\boldsymbol{a}) + N_{S}v_{int_S}(\boldsymbol{bS},\boldsymbol{a})\right)$

#### Transcription
$$
\begin{flalign*}
&\omega_q(\boldsymbol{a}) = Q \cdot w_{q}\frac{\boldsymbol{a}}{\theta_{q}+\boldsymbol{a}} && \\
&\omega_L(\boldsymbol{L},\boldsymbol{a}) = \boldsymbol{L} \cdot w_{L}\frac{\boldsymbol{a}}{\theta_{L}+\boldsymbol{a}} && \\
&\omega_S(\boldsymbol{S},\boldsymbol{a}) = \boldsymbol{S} \cdot w_{S}\frac{\boldsymbol{a}}{\theta_{S}+\boldsymbol{a}}
\end{flalign*}
$$

#### Translation
$$
\begin{flalign*}
&v_q(\boldsymbol{cq},\boldsymbol{a}) = \frac{\gamma_{\max_q}}{N_{q}}\boldsymbol{cq}\frac{\boldsymbol{a}}{K_{\gamma_q}+\boldsymbol{a}} && \\
&v_L(\boldsymbol{cL},\boldsymbol{a}) = \frac{\gamma_{\max_L}}{N_{L}/3}\boldsymbol{cL}\frac{\boldsymbol{a}}{K_{\gamma_L}+\boldsymbol{a}}
\end{flalign*}
$$

#### Integration
$$
\begin{flalign*}
&v_{int_L}(\boldsymbol{bL},\boldsymbol{01},\boldsymbol{a}) = \frac{\chi_{\max_L}}{N_{L}} \boldsymbol{bL} \frac{\boldsymbol{a}}{K_{\chi_L}+\boldsymbol{a}} \cdot \frac{K_{L}\boldsymbol{O_1}}{1+K_{L}\boldsymbol{O_1}} && \\
&v_{int_S}(\boldsymbol{bS},\boldsymbol{a}) = \frac{\chi_{\max_S}}{N_{S}} \boldsymbol{bS} \frac{\boldsymbol{a}}{K_{\chi_S}+\boldsymbol{a}}
\end{flalign*}
$$

#### Free ribosomes
$\quad f_{rib}(\boldsymbol{cq},\boldsymbol{cL}) = r_{tot} - \boldsymbol{cq} - \boldsymbol{cL}$

### Table with model parameters (reference model)
| Designation       | Description                                                               | Value             | Unit        | Source  
| :---              | :---                                                                      | :---              | :---        | :---      
| $a_{wt}$          | characteristic number of ATP molecules in HeLa cell                       | $5.33·10^9$       | m.p.c.      | $BNID:104449$
| $L_{wt}$          | characteristic number of active LINE-1 in HeLa cell                       | 1064              | m.p.c.      | $UCSS$
| $S_{wt}$          | characteristic number of active SINE in HeLa cell                         | 13243             | m.p.c.      | $UCSS$
| $N_g$             | total genome size                                                         | $3.08·10^9$       | bp          | $BNID:101484$
| $Q$               | number of housekeeping genes (genes of $\boldsymbol{q}$)                  | 3804              |             | $[3]$
| $N_q$             | median HeLa cell protein length ($\boldsymbol{q}$)                        | 431               | aa          | $[4]$
| $N_Q$             | median HeLa gen length (genes of $\boldsymbol{q}$)                        | 1300              | bp          | $[4]^2$
| $N_L$             | LINE-1 length                                                             | 6000              | bp          | $[5]$
| $N_S$             | SINE length                                                               | 300               | bp          | $[6]$
| $N_{aa}$          | number of ATP molecules for adding one aa                                 | 5                 | m.p.c.      | $[7]$
| $N_{nt}$          | number of ATP molecules for adding one nucleotide                         | 15                | m.p.c.      | $[7]$
| $\tau$            | HeLa cell cycle duration                                                  | 1320              | min         | $BNID:109393$
| $r_{tot}$         | total number of ribosomes                                                 | $9.5·10^6$        | m.p.c.      | $BNID:107347$
| $A_0$             | energy supply                                                             | $1.31·10^{10}$    | m.p.c./min  | $BNID:110879^1$
| $\chi_{\max_L}$   | maximal integration elongation rate of one LINE-1 transposon base pair    | 840               | bp/min      | $[8]^3$
| $\chi_{\max_S}$   | maximal integration elongation rate of one SINE transposon base pair      | 840               | bp/min      | $[8]$
| $K_{\chi_L}$      | integration elongation threshold of one LINE-1 transposon                 | $1.1·10^7$        | m.p.c.      | $[8]^4$
| $K_{\chi_S}$      | integration elongation threshold of one SINE transposon                   | $1.1·10^7$        | m.p.c.      | $[8]$
| $K_L$             | association constant of LINE-1 mRNA with ORF1p                            | $2.24·10^{-3}$    | 1/m.p.c.    | $[9]^5$
| $w_q$             | maximal transcription rate of one housekeeping gene                       | 4.64              | m.p.c./min  | $BNID:111721^6$
| $w_L$             | maximal transcription rate of one LINE-1                                  | 1                 | m.p.c./min  | $BNID:111721$
| $w_S$             | maximal transcription rate of one SINE                                    | 20                | m.p.c./min  | $BNID:111721$
| $\theta_q$        | transcription threshold of one housekeeping gene                          | $3.8·10^9$        | m.p.c.      | $BNID:111027^7$
| $\theta_L$        | transcription threshold of one LINE-1                                     | $3.8·10^9$        | m.p.c.      | $BNID:111027$
| $\theta_S$        | transcription threshold of one SINE                                       | $3.8·10^9$        | m.p.c.      | $BNID:111027$
| $\gamma_{\max_q}$ | maximal translation rate of one aa from q-RNA                             | 300               | aa/min      | $BNID:104598^8$
| $\gamma_{\max_L}$ | maximal translation rate of one aa from LINE-1 RNA                        | 300               | aa/min      | $BNID:104598$
| $K_{\gamma_q}$    | translation threshold of one $\boldsymbol{q}$-RNA                         | 25900             | m.p.c.      | $[10]^9$
| $K_{\gamma_L}$    | translation threshold of one LINE-1 RNA                                   | 25900             | m.p.c.      | $[10]$
| $k_{bq}$          | $\boldsymbol{cq}$-complexes (RNA+ribosome) binding rate constant          | $5·10^{-8}$       | 1/min       | $*$
| $k_{bL}$          | $\boldsymbol{cL}$-complexes (RNA+ribosome) binding rate constant          | $5·10^{-8}$       | 1/min       | $*$
| $k_{uq}$          | $\boldsymbol{cq}$-complexes unbinding rate constant                       | 0.01              | 1/min       | $*$
| $k_{uL}$          | $\boldsymbol{cL}$-complexes unbinding rate constant                       | 0.01              | 1/min       | $*$
| $k_{sub_S}$       | substitution of mL on mS rate constant in complex with ORF2p              | $5·10^{-8}$       | 1/min       | $*$
| $k_{sub_L}$       | substitution of mS on mL rate constant in complex with ORF2p              | $5·10^{-6}$       | 1/min       | $*$
| $d_{mq}$          | $\boldsymbol{q}$-RNAs degradation rate constant                           | $1.15·10^{-3}$    | 1/min       | $BNID:104747$
| $d_{cq}$          | $\boldsymbol{cq}$-complexes degradation rate constant                     | $1.55·10^{-3}$    | 1/min       | $*$
| $d_{q}$           | $\boldsymbol{q}$ proteins degradation rate constant                       | $5.67·10^{-4}$    | 1/min       | $BNID:112253$
| $d_{mL}$          | LINE-1 RNAs degradation rate constant                                     | $1.15·10^{-3}$    | 1/min       | $BNID:104747$
| $d_{mS}$          | SINE RNAs degradation rate constant                                       | $1.15·10^{-3}$    | 1/min       | $BNID:104747$
| $d_{cL}$          | $\boldsymbol{cL}$-complexes degradation rate constant                     | $1.55·10^{-3}$    | 1/min       | $*$
| $d_{O1}$          | ORF1p degradation rate constant                                           | $5.67·10^{-4}$    | 1/min       | $BNID:112253$
| $d_{bL}$          | $\boldsymbol{bL}$-complexes degradation rate constant                     | $5.67·10^{-4}$    | 1/min       | $BNID:112253$
| $d_{bS}$          | $\boldsymbol{bS}$-complexes degradation rate constant                     | $5.67·10^{-4}$    | 1/min       | $BNID:112253$
| $\lambda_{a}$     | ATP molecules degradation rate constant                                   | 1.47              | 1/min       | $*$
| $\lambda_{L}$     | LINE-1 deactivation rate constant                                         | 0.37              | 1/min       | $*$
| $\lambda_{S}$     | SINE deactivation rate constant                                           | 1.18              | 1/min       | $*$
| $V_{cell}$        | cell volume                                                               | 3700              | $μm^3$      | $BNID:105879$

### Table with modified parameters (relative to the reference model)
| Designation       | Description                                                               | Value             | Unit        | Source  
| :---              | :---                                                                      | :---              | :---        | :---      
| $\lambda_{a}$     | ATP molecules degradation rate constant                                   | $2.04·10^{-3}$    | 1/min       | $*$
| $\lambda_{L}$     | LINE-1 deactivation rate constant                                         | 0.53              | 1/min       | $*$
| $\lambda_{S}$     | SINE deactivation rate constant                                           | 2.47              | 1/min       | $*$

### Stationary solutions

#### Reference model (ATP value from literature)
| Designation       | Active transposons    | Knock out SINE        | Knock out SINE & LINE-1   | Active transposons, $\lambda_{L}/2$   | Active transposons, $\lambda_{L}/2$ & $\lambda_{S}/2$ |
| :---              | :---                  | :---                  | :---                      | :---                                  | :---                                                  |
| $\boldsymbol{a}$  | 5.32028879e+09        | 2.29529555e+09        | 5.96205458e+09            | 3.60500385e+09                        | 4.27240077e+09                                        |
| $\boldsymbol{mq}$ | 5.43892043e+06        | 3.36327244e+06        | 5.69819032e+06            | 4.52734948e+06                        | 4.89876923e+06                                        |
| $\boldsymbol{cq}$ | 2.60750605e+06        | 1.79282544e+06        | 2.72706015e+06            | 2.18479470e+06                        | 2.39235807e+06                                        |
| $\boldsymbol{q}$  | 3.20098867e+09        | 2.20086818e+09        | 3.34775568e+09            | 2.68205969e+09                        | 2.93686833e+09                                        |
| $\boldsymbol{mL}$ | 5.13241818e+04        | 6.99550598e+04        | 0                         | 2.29815368e+05                        | 9.18023587e+04                                        |
| $\boldsymbol{mS}$ | 1.21841508e+08        | 0                     | 0                         | 4.55392447e+08                        | 4.30559873e+08                                        |
| $\boldsymbol{cL}$ | 1.07775213e+05        | 1.63334967e+05        | 0                         | 4.85769072e+05                        | 1.96370953e+05                                        |
| $\boldsymbol{O_1}$| 2.78154330e+07        | 1.75569280e+05        | 0                         | 1.27407054e+08                        | 5.15602878e+07                                        |
| $\boldsymbol{bL}$ | 2.82625278e+03        | 1.75569280e+05        | 0                         | 4.47819184e+03                        | 1.58108214e+03                                        |
| $\boldsymbol{bS}$ | 5.64255248e+03        | 0                     | 0                         | 2.58724875e+04                        | 1.04653946e+04                                        |
| $\boldsymbol{L}$  | 1.06716935e+03        | 6.59470834e+04        | 0                         | 3.37858095e+03                        | 1.19341160e+03                                        |
| $\boldsymbol{S}$  | 1.33614820e+04        | 0                     | 0                         | 6.12055856e+04                        | 4.95387336e+04                                        |

#### Nota bene
Stationary solutions and changable parameters for other models can be found in **.py** scripts in **stat_sols/** folder. 


### Petersen Matrix (transposed stoichiometric matrix)
| Process                                           | a                                 | mq   | cq   | q    | mL   | mS   | cL   | O1   | bL   | bS   | L    | S    |
| :---                                              | :---                              | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Energy supply                                     | 1                                 |      |      |      |      |      |      |      |      |      |      |      |
| Replication                                       | -1                                |      |      |      |      |      |      |      |      |      |      |      |
| Transcription of housekeeping genes               | -1 $N_{nt}N_{Q}$                  | 1    |      |      |      |      |      |      |      |      |      |      |
| Transcription of LINE-1                           | -1 $N_{nt}N_{L}$                  |      |      |      | 1    |      |      |      |      |      |      |      |
| Transcription of SINE                             | -1 $N_{nt}N_{S}$                  |      |      |      |      | 1    |      |      |      |      |      |      |
| Ribosome binding with housekeeping genes RNA      |                                   | -1   | 1    |      |      |      |      |      |      |      |      |      |
| Ribosome binding with LINE-1 RNA                  |                                   |      |      |      | -1   |      | 1    |      |      |      |      |      |
| Ribosome unbinding housekeeping genes RNA         |                                   | 1    | -1   |      |      |      |      |      |      |      |      |      |
| Ribosome unbinding LINE-1 RNA                     |                                   |      |      |      | 1    |      | -1   |      |      |      |      |      |
| Translation of housekeeping genes RNAs            | -1 $N_{aa}N_{q}$                  | 1    | -1   | 1    |      |      |      |      |      |      |      |      |
| Translation of LINE-1 RNAs                        | -1 $N_{aa}\frac{N_{L}}{3}\quad$   |      |      |      |      |      | -1   | 1    | 1    |      |      |      |
| Replacement $mL$ --> $mS$ in complex with $O_2$   |                                   |      |      |      | 1    | -1   |      |      | -1   | 1    |      |      |
| Replacement $mS$ --> $mL$ in complex with $O_2$   |                                   |      |      |      | -1   | 1    |      |      | 1    | -1   |      |      |
| Integration of LINE-1                             | -1 $N_{nt}N_{L}$                  |      |      |      |      |      |      | -1   | -1   |      | 1    |      |
| Integration of SINE                               | -1 $N_{nt}N_{S}$                  |      |      |      |      |      |      |      |      | -1   |      | 1    |
| ATP degradation                                   | -1                                |      |      |      |      |      |      |      |      |      |      |      |
| $mq$ degradation                                  |                                   | -1   |      |      |      |      |      |      |      |      |      |      |
| $cq$ degradation                                  |                                   |      | -1   |      |      |      |      |      |      |      |      |      |
| housekeeping genes degradation                    |                                   |      |      | -1   |      |      |      |      |      |      |      |      |
| $mL$ degradation                                  |                                   |      |      |      | -1   |      |      |      |      |      |      |      |
| $mS$ degradation                                  |                                   |      |      |      |      | -1   |      |      |      |      |      |      |
| $cL$-complex degradation                          |                                   |      |      |      |      |      | -1   |      |      |      |      |      |
| $O_1$ degradation                                 |                                   |      |      |      |      |      |      | -1   |      |      |      |      |
| $bL$-complex degradation                          |                                   |      |      |      |      |      |      |      | -1   |      |      |      |
| $bS$-complex degradation                          |                                   |      |      |      |      |      |      |      |      | -1   |      |      |
| LINE-1 degradation                                |                                   |      |      |      |      |      |      |      |      |      | -1   |      |
| SINE degradation                                  |                                   |      |      |      |      |      |      |      |      |      |      | -1   |
