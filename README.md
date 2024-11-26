# TE_model
This repository contains code with a model describing the basic biochemical processes in a cell in conditions of limited bioenergetic resources and competition for ATP and ribosomes from retrotransposons.
The model separately considers the processes associated with housekeeping genes and two types of non-LTR retrotransposons: LINEs and SINEs.

## Languages
We use Jupiter notebooks for interactive charting.

## Description of the mathematical model

### Dynamical variables:
a	- amount of available energy (ATP)
mq	- the RNA of housekeeping genes
cq	- complexes of mq with ribosomes
q	- proteins translated from cq
mL	- RNA of LINE-1
mS	- RNA of SINE
cL	- complexes of mL with ribosomes
O1	- ORF1p proteins translated from cL (needed for new LINE-1 integration)
bL	- complexes of ORF2p with mL
bS	- complexes of ORF2p proteins with mS
L	- LINE-1 transposons
S	- SINE transposons

### Ordinary Differential Equations
$\begin{aligned}
\frac{da}{dt}=A_{0}-\lambda_{a}& a- v_{repl}(\boldsymbol{a},\boldsymbol{b}\boldsymbol{L},\boldsymbol{O}\boldsymbol{1},\boldsymbol{b}\boldsymbol{S})- N_{nt}\left(N_{Q}\omega_{q}(\boldsymbol{a})+N_{S}\omega_{S}(\boldsymbol{S},\boldsymbol{a})+N_{L}\omega_{L}(\boldsymbol{L},\boldsymbol{a})\right)- N_{aa}\left(N_{q}v_{q}(\boldsymbol{c}\boldsymbol{q},\boldsymbol{a})+\right) \\
&\frac{N_L}{3}v_L(\boldsymbol{cL},\boldsymbol{a})\Big)- N_{nt}(N_Lv_{\mathrm{int}L}(\boldsymbol{bL},\boldsymbol{O1},\boldsymbol{a})+N_Sv_{\mathrm{int}S}(\boldsymbol{bS},\boldsymbol{a})\Big)& \left(1\right) \\
&\frac{dmq}{dt}=\omega_{q}(\boldsymbol{a})-k_{bq}f_{rib}(\boldsymbol{cq},\boldsymbol{cL})\boldsymbol{mq}+v_{q}(\boldsymbol{cq},\boldsymbol{a})+k_{uq}\boldsymbol{cq}-d_{mq}\boldsymbol{mq}& (2) \\
&\frac{dcq}{dt}=k_{bq}f_{rib}\left(cq,cL\right)mq-k_{uq}cq-v_{q}\left(cq,a\right)-d_{cq}cq& \left(3\right) \\
&\frac{dq}{dt}=v_{q}(cq,a)-d_{q}q& (4) \\
&\frac{dmL}{dt}=\omega_{L}(L,a)-k_{bL}f_{rib}(cq,cL)mL+k_{uL}cL+k_{\mathrm{sub}S}mS bL-k_{\mathrm{sub}L}mL bS-d_{mL}mL& \left(5\right) \\
&\frac{dmS}{dt}=\omega_{S}(S,a)-k_{\mathrm{sub}S}mS bL+k_{\mathrm{sub}L}mL bS-d_{mS}mS& (6) \\
&\frac{dcL}{dt}=k_{bL}f_{rib}(cq,cL)mL-k_{uL}cL-v_{L}(cL,a)-d_{cL}cL& \left(7\right) \\
&\frac{dO1}{dt}=v_{L}(cL,a)-v_{\mathrm{int}L}(bL,01,a)-d_{O1}O1& \left(8\right) \\
&\frac{dbL}{dt}=v_{L}(cL,a)-v_{\mathrm{int}L}(bL,O1,a)-k_{\mathrm{sub}S}\boldsymbol{mS}\boldsymbol{bL}+k_{\mathrm{sub}L}\boldsymbol{mL}\boldsymbol{bS}-d_{bL}\boldsymbol{bL}& \left(9\right) \\
&\frac{dbS}{dt}=k_{\mathrm{sub}S}\boldsymbol{mSbL}-k_{\mathrm{sub}L}\boldsymbol{mL}\boldsymbol{bS}-v_{\mathrm{int}S}(\boldsymbol{bS},\boldsymbol{a})-d_{bS}\boldsymbol{bS}& (10) \\
&\frac{dL}{dt}=v_{\mathrm{int}L}(bL,01,a)-\lambda_{L}L& (11) \\
&\frac{dS}{dt}=v_{\mathrm{int}S}(bS,a)-\lambda_{S}S& (12) 
\end{aligned}$