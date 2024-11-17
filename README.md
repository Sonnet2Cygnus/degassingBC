# Degassing Boundary Condition for OpenFOAM10

The main reference equations are as follows:

Mass source:
$$
S_D=-\rho_\text{cell}\cdot\text{VOF}_\text{cell}/\Delta t
$$

The momentum source for the continuous phase：
$$
S_{RP}=-\rho_{\mathrm{sub}}\cdot\mathrm{VOF_{sub}}/\Delta t\cdot u_{\mathrm{prim}}
$$

The momentum source for the dispersed phase:
$$
S_{RD}=-\rho_{\mathrm{sub}}\cdot\mathrm{VOF_{sub}}/\Delta t\cdot u_{\mathrm{sub}}
$$
