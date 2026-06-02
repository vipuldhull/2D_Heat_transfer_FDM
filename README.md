# 2D_Heat_transfer_FDM
Solver based in python using finite difference method to solve 2D steady state heat transfer
# Stage 1 — Derivation of Finite Difference Scheme

Derive a **second-order accurate Finite Difference (FD)** scheme for the **2D steady-state heat equation** with **constant thermal conductivity** $\lambda$.  
Assume **Dirichlet boundary conditions** are applied on all sides of the domain.

1. Start from the continuous governing equation  
3. Apply **central differencing** to discretize spatial derivatives.  
## Solution — Stage 1

For constant thermal conductivity $\lambda$, the steady heat equation becomes

$$
\lambda \left(\frac{\partial^2 T}{\partial x^2}+\frac{\partial^2 T}{\partial y^2}\right)=0
$$

Since $\lambda$ is constant and not zero, we can divide by $\lambda$:

$$
\frac{\partial^2 T}{\partial x^2}+\frac{\partial^2 T}{\partial y^2}=0
$$

Using second-order central differences at an inner node $(i,j)$:

$$
\frac{\partial^2 T}{\partial x^2}\Bigg|_{i,j}
\approx
\frac{T_{i-1,j}-2T_{i,j}+T_{i+1,j}}{\Delta x^2}
$$

$$
\frac{\partial^2 T}{\partial y^2}\Bigg|_{i,j}
\approx
\frac{T_{i,j-1}-2T_{i,j}+T_{i,j+1}}{\Delta y^2}
$$

Therefore the finite difference equation is

$$
\frac{T_{i-1,j}-2T_{i,j}+T_{i+1,j}}{\Delta x^2}
+
\frac{T_{i,j-1}-2T_{i,j}+T_{i,j+1}}{\Delta y^2}
=0
$$

For a square grid where $\Delta x=\Delta y$, this simplifies to

$$
T_{i+1,j}+T_{i-1,j}+T_{i,j+1}+T_{i,j-1}-4T_{i,j}=0
$$

or

$$
T_{i,j}=\frac{T_{i+1,j}+T_{i-1,j}+T_{i,j+1}+T_{i,j-1}}{4}
$$

So the temperature at an inner node is the average of its four neighbours.
