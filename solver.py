import numpy as np
import matplotlib.pyplot as plt

class Source:
    def __init__(self, location, is_active, q_dot):
        self.location  = location   # [x, y] physical coordinates
        self.is_active = is_active  # bool
        self.q_dot     = q_dot      # intensity in [W/m]


class SteadyHeat2D:
    def __init__(self, Lx, Ly, dimX, dimY, source):
        self.l    = Lx #domain length in x
        self.h    = Ly #domain length in y
        self.dimX = dimX #number of nodes in x
        self.dimY = dimY #number of nodes in y

        self.source_location = source.location
        self.source_active   = source.is_active
        self.source_q        = source.q_dot

        # Grid spacing
        self.dx = Lx / (dimX - 1)
        self.dy = Ly / (dimY - 1)

        # Linear system A * T_flat = b
        N    = dimX * dimY
        self.A = np.zeros((N, N))
        self.b = np.zeros(N)
        self.set_inner()

    def _idx(self, i, j):
        """Convert 2D grid index (i, j) to 1D matrix row index."""
        return i * self.dimY + j

    # build linear system
    def set_inner(self):
        """
        Assemble coefficient matrix A and RHS vector b

        If a point source is active, add its contribution
        to the corresponding node in the RHS vector.
        """
        dx, dy = self.dx, self.dy

        for i in range(1, self.dimX - 1):
            for j in range(1, self.dimY - 1):
                k = self._idx(i, j)
                # Center coefficient
                self.A[k, k]                   = -2/dx**2 - 2/dy**2
                # West neighbor
                self.A[k, self._idx(i-1, j)]   =  1/dx**2
                # East neighbor
                self.A[k, self._idx(i+1, j)]   =  1/dx**2
                # South neighbor
                self.A[k, self._idx(i,   j-1)] =  1/dy**2
                # North neighbor
                self.A[k, self._idx(i,   j+1)] =  1/dy**2
                self.b[k] = 0.0

        # --- Point source ---
        if self.source_active:
            x_s, y_s = self.source_location
            i_s = int(round(x_s / dx))
            j_s = int(round(y_s / dy))
            i_s = max(1, min(self.dimX - 2, i_s))  # clamp to interior
            j_s = max(1, min(self.dimY - 2, j_s))
            k_s = self._idx(i_s, j_s)
            self.b[k_s] -= self.source_q / (dx * dy)


    # Set boundary condition )
    # south 
    def set_south(self, bc_type, T_d=0.0, q=0.0, alpha=0.0, T_inf=0.0):
        """
        Apply south boundary condition.

        Parameters
        ----------
        bc_type : str
            Type of boundary condition ('D', 'N', 'R')
        T_d : float
            Prescribed temperature (for Dirichlet)
        q : float
            Prescribed heat flux (for Neumann)
        alpha : float
            Convection coefficient (for Robin)
        T_inf : float
            Ambient temperature (for Robin)
        """
        dx, dy = self.dx, self.dy
        for i in range(1, self.dimX - 1):
            k = self._idx(i, 0)
            self.A[k, :] = 0.0;  self.b[k] = 0.0

            if bc_type == 'D':
                self.A[k, k] = 1.0
                self.b[k]    = T_d

            elif bc_type == 'N':
                self.A[k, k]                 = -2/dx**2 - 2/dy**2
                self.A[k, self._idx(i-1, 0)] =  1/dx**2
                self.A[k, self._idx(i+1, 0)] =  1/dx**2
                self.A[k, self._idx(i,   1)] =  2/dy**2
                self.b[k] = 2*q/dy

            elif bc_type == 'R':
                self.A[k, k]                 = -2/dx**2 - 2/dy**2 - 2*alpha/dy
                self.A[k, self._idx(i-1, 0)] =  1/dx**2
                self.A[k, self._idx(i+1, 0)] =  1/dx**2
                self.A[k, self._idx(i,   1)] =  2/dy**2
                self.b[k] = -2*alpha*T_inf/dy

    # North boundary 
    def set_north(self, bc_type, T_d=0.0, q=0.0, alpha=0.0, T_inf=0.0):
        """
        Apply north boundary condition.
        Same parameters as set_south.
        """
        dx, dy = self.dx, self.dy
        jN = self.dimY - 1
        for i in range(1, self.dimX - 1):
            k = self._idx(i, jN)
            self.A[k, :] = 0.0;  self.b[k] = 0.0

            if bc_type == 'D':
                self.A[k, k] = 1.0
                self.b[k]    = T_d

            elif bc_type == 'N':
                self.A[k, k]                    = -2/dx**2 - 2/dy**2
                self.A[k, self._idx(i-1, jN)]   =  1/dx**2
                self.A[k, self._idx(i+1, jN)]   =  1/dx**2
                self.A[k, self._idx(i,   jN-1)] =  2/dy**2
                self.b[k] = 2*q/dy

            elif bc_type == 'R':
                self.A[k, k]                    = -2/dx**2 - 2/dy**2 - 2*alpha/dy
                self.A[k, self._idx(i-1, jN)]   =  1/dx**2
                self.A[k, self._idx(i+1, jN)]   =  1/dx**2
                self.A[k, self._idx(i,   jN-1)] =  2/dy**2
                self.b[k] = -2*alpha*T_inf/dy

    # West 
    def _west_row(self, j, bc_type, T_d, q, alpha, T_inf):
        """Fill one row for the west boundary at column j."""
        dx, dy = self.dx, self.dy
        k = self._idx(0, j)
        self.A[k, :] = 0.0;  self.b[k] = 0.0

        if bc_type == 'D':
            self.A[k, k] = 1.0;  self.b[k] = T_d;  return

        # --- x-part 
        if bc_type == 'N':
            self.A[k, k]              = -2/dx**2
            self.A[k, self._idx(1,j)] =  2/dx**2
            self.b[k] = 2*q/dx
        elif bc_type == 'R':
            self.A[k, k]              = -2/dx**2 - 2*alpha/dx
            self.A[k, self._idx(1,j)] =  2/dx**2
            self.b[k] = -2*alpha*T_inf/dx

        # --- y-part ---  corner zero flux approximation
        if 0 < j < self.dimY - 1:          # interior j: both neighbors present
            self.A[k, k]                    -= 2/dy**2
            self.A[k, self._idx(0, j-1)]     = 1/dy**2
            self.A[k, self._idx(0, j+1)]     = 1/dy**2
        elif j == 0:                        # SW corner: T_{0,-1} = T_{0,0}
            self.A[k, k]                    -= 1/dy**2
            self.A[k, self._idx(0, 1)]       = 1/dy**2
        else:                               # NW corner: T_{0,dimY} = T_{0,dimY-1}
            self.A[k, k]                    -= 1/dy**2
            self.A[k, self._idx(0, self.dimY-2)] = 1/dy**2

    # East 
    def _east_row(self, j, bc_type, T_d, q, alpha, T_inf):
        """Fill one row for the east boundary at column j."""
        dx, dy = self.dx, self.dy
        iE = self.dimX - 1
        k  = self._idx(iE, j)
        self.A[k, :] = 0.0;  self.b[k] = 0.0

        if bc_type == 'D':
            self.A[k, k] = 1.0;  self.b[k] = T_d;  return

        if bc_type == 'N':
            self.A[k, k]                  = -2/dx**2
            self.A[k, self._idx(iE-1, j)] =  2/dx**2
            self.b[k] = 2*q/dx
        elif bc_type == 'R':
            self.A[k, k]                  = -2/dx**2 - 2*alpha/dx
            self.A[k, self._idx(iE-1, j)] =  2/dx**2
            self.b[k] = -2*alpha*T_inf/dx

        # --- y-part ---for corner zero flux approximation 
        if 0 < j < self.dimY - 1:
            self.A[k, k]                     -= 2/dy**2
            self.A[k, self._idx(iE, j-1)]     = 1/dy**2
            self.A[k, self._idx(iE, j+1)]     = 1/dy**2
        elif j == 0:
            self.A[k, k]                     -= 1/dy**2
            self.A[k, self._idx(iE, 1)]       = 1/dy**2
        else:
            self.A[k, k]                     -= 1/dy**2
            self.A[k, self._idx(iE, self.dimY-2)] = 1/dy**2

    def set_west(self, bc_type, T_d=0.0, q=0.0, alpha=0.0, T_inf=0.0):
        """Apply west BC (i=0) for all j including SW and NW corners."""
        for j in range(self.dimY):
            self._west_row(j, bc_type, T_d, q, alpha, T_inf)

    def set_east(self, bc_type, T_d=0.0, q=0.0, alpha=0.0, T_inf=0.0):
        """Apply east BC (i=dimX-1) for all j including SE and NE corners."""
        for j in range(self.dimY):
            self._east_row(j, bc_type, T_d, q, alpha, T_inf)

    # Solve
    def solve(self):
        """
        Solve the assembled linear system A·T = b
        using np.lingalg.solve()
        and reshape the temperature field to 2D.
        """
        T_flat = np.linalg.solve(self.A, self.b)
        return T_flat.reshape(self.dimX, self.dimY)

    # Visualise
    # ------------------------------------------------------------------
    def plot(self, T, title='2D Steady-State Heat Conduction'):
        """
        Visualize the temperature field as a 2D contour / heatmap.
        """

        x = np.linspace(0, self.l, self.dimX)
        y = np.linspace(0, self.h, self.dimY)
        X, Y = np.meshgrid(x, y)

        plt.figure(figsize=(8, 6))

        # --- Heatmap only ---
        cp = plt.contourf(X, Y, T.T, levels=50, cmap='hot')
        plt.colorbar(cp, label='Temperature [K]')

        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.title(title)

        plt.tight_layout()
        plt.show()

        print(f'T_min = {T.min():.2f} K   |   T_max = {T.max():.2f} K')
