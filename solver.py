import numpy as np
import matplotlib.pyplot as plt

class source:
  def _init_(self, location, is_active, q_dot):
    self.location = location
    self.is_active = is_active
    self.q_dot = q_dot

class Steadyheat2D:
  def _init_(self, Lx, Ly, DimX, DimY, source):
    self.L = Lx
    self.H = Ly
    self.dimX = DimX
    self.dimY = DimY
    self.source = source

    
