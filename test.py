# ================================================================
# Test Case 1 — All Dirichlet Boundary Conditions
# ================================================================

# Domain and grid setup
Lx = 1.0
Ly = 1.0
dimX = 101
dimY = 101

source = Source(location=[0.5, 0.5], is_active=False, q_dot=0.0)

heat = SteadyHeat2D(Lx, Ly, dimX, dimY, source)

heat.set_south('D', T_d=300)   # Bottom — hot
heat.set_north('N', T_d=100)   # Top
heat.set_east( 'D', T_d=100)   # Right
heat.set_west( 'D', T_d=100)   # Left

T = heat.solve()               # shape = (dimX, dimY)
heat.plot(T, 'TC1: All Dirichlet BCs  —  T_south=300 K, others=100 K')

# ================================================================
# Test Case 2 — Mixed Boundary Conditions - Original Case
# ================================================================

# Domain and grid setup
Lx = 1.0
Ly = 1.0
dimX = 101
dimY = 101

source = Source(location=[0.5, 0.5], is_active=False, q_dot=0.0)

heat = SteadyHeat2D(Lx, Ly, dimX, dimY, source)
heat.set_south('N', q=0.0)                      # Insulated (zero flux)
heat.set_north('R', alpha=1.0, T_inf=100)       # Convective
heat.set_east( 'R', alpha=1.0, T_inf=100)       # Convective
heat.set_west( 'D', T_d=300)                    # Fixed temperature

T = heat.solve()
heat.plot(T, 'TC2: Dirichlet West | Neumann South (q=0) | Robin North & East')

# ================================================================
# Test Case 3 — Mixed Boundary Conditions with Point Source Active
# ================================================================

Lx, Ly = 1.0, 1.0
dimX = 101
dimY = 101

# Point source at (x=0.5, y=0.95), near the north Robin boundary
source = Source(location=[0.5, 0.95], is_active=True, q_dot=200)

heat = SteadyHeat2D(Lx, Ly, dimX, dimY, source)
heat.set_south('N', q=0.0)
heat.set_north('R', alpha=1.0, T_inf=100)
heat.set_east( 'R', alpha=1.0, T_inf=100)
heat.set_west( 'D', T_d=300)

T = heat.solve()
heat.plot(T, 'TC3: Mixed BCs + Point Source at (0.5, 0.95)  [q_dot = 200 W/m]')
