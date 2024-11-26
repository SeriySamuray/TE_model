# TE_model
This repository contains code with a model describing the basic biochemical processes in a cell in conditions of limited bioenergetic resources and competition for ATP and ribosomes from retrotransposons.
The model separately considers the processes associated with housekeeping genes and two types of non-LTR retrotransposons: LINEs and SINEs.

## Languages
We use Jupiter notebooks for interactive charting.

## Description of the mathematical model

### Dynamical variables:
```math
\begin{flalign*}
\boldsymbol{a}  &\text{- amount of available energy (ATP)}&&\\
\boldsymbol{mq} &\text{- the RNA of housekeeping genes}&&\\
\boldsymbol{cq} &\text{- complexes of \textbf{mq} with ribosomes}&&\\
\boldsymbol{q}  &\text{- proteins translated from \textbf{cq}}&&\\
\boldsymbol{mL} &\text{- RNA of LINE-1}&&\\
\boldsymbol{mS} &\text{- RNA of SINE}&&\\
\boldsymbol{cL} &\text{- complexes of \textbf{mL} with ribosomes}&&\\
\boldsymbol{O_1}&\text{- ORF1p proteins translated from \textbf{cL} (needed for new LINE-1 integration)}&&\\
\boldsymbol{bL} &\text{- complexes of ORF2p with \textbf{mL}}&&\\
\boldsymbol{bS} &\text{- complexes of ORF2p proteins with \textbf{mS}}&&\\
\boldsymbol{L}  &\text{- LINE-1 transposons}&&\\
\boldsymbol{S}  &\text{- SINE transposons}&&
\end{flalign*}
```

### Ordinary Differential Equations
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

& 4. \quad \frac{\boldsymbol{dq}}{\boldsymbol{dt}} = v_{q}(\boldsymbol{cq},\boldsymbol{a}) - d_{q}\boldsymbol{q} && \\

& 5. \quad \frac{\boldsymbol{dmL}}{\boldsymbol{dt}} = \omega_{L}(\boldsymbol{L},\boldsymbol{a}) - k_{bL}f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\:\boldsymbol{mL} +
   k_{uL}\boldsymbol{cL} + k_{\mathrm{sub}S}mS bL - k_{sub_L}\boldsymbol{mL}\boldsymbol{bS} - d_{mL}\boldsymbol{mL} && \\

& 6. \quad \frac{\boldsymbol{dmS}}{\boldsymbol{dt}} = \omega_{S}(\boldsymbol{S},\boldsymbol{a}) - k_{sub_S}\boldsymbol{mS}\boldsymbol{bL} + k_{sub_L}\boldsymbol{mL}\boldsymbol{bS} - d_{mS}\boldsymbol{mS} && \\

& 7. \quad \frac{\boldsymbol{dcL}}{\boldsymbol{dt}} = k_{bL}f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\boldsymbol{mL} - k_{uL}\boldsymbol{cL} - v_{L}(\boldsymbol{cL},\boldsymbol{a}) - d_{cL}\boldsymbol{cL} && \\

& 8. \quad \frac{\boldsymbol{dO1}}{\boldsymbol{dt}} = v_{L}(cL,a) - v_{int_L}(\boldsymbol{bL},\boldsymbol{01},\boldsymbol{a}) - d_{O1}\boldsymbol{O1} && \\

& 9. \quad \frac{\boldsymbol{dbL}}{\boldsymbol{dt}} = v_{L}(cL,a) - v_{int_L}(\boldsymbol{bL},\boldsymbol{O1},\boldsymbol{a}) - k_{sub_S}\boldsymbol{mS}\boldsymbol{bL} + k_{sub_L}\boldsymbol{mL}\boldsymbol{bS} - d_{bL}\boldsymbol{bL} && \\

& 10.\quad \frac{\boldsymbol{dbS}}{\boldsymbol{dt}} = k_{sub_S}\boldsymbol{mS}\boldsymbol{bL} - k_{sub_L}\boldsymbol{mL}\boldsymbol{bS} - v_{int_S}(\boldsymbol{bS},\boldsymbol{a}) - d_{bS}\boldsymbol{bS} && \\

& 11.\quad \frac{\boldsymbol{dL}}{\boldsymbol{dt}} = v_{int_L}(\boldsymbol{bL},\boldsymbol{01},\boldsymbol{a}) - \lambda_{L}\boldsymbol{L} && \\

& 12.\quad \frac{\boldsymbol{dS}}{\boldsymbol{dt}} = v_{int_S}(\boldsymbol{bS},\boldsymbol{a}) - \lambda_{S}\boldsymbol{S} && 
\end{flalign*}
```


