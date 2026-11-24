import numpy as np
import pandas as pd

def sim(para):
    
  # Unpack parameters
  nubar = para['nubar']
  theta = para['theta']
  delta = para['delta']
  boundary = para['boundary']
  ngen = para['ngen']
  ninit = para['ninit']

  # Initialize variables
  parent = np.full(boundary, np.nan)
  nup = np.full(boundary, np.nan)
  nu = np.full(boundary, np.nan)
  onp = np.full(boundary, np.nan)
  on = np.full(boundary, np.nan)
  gen = np.full(boundary, np.nan)
  
  boundary_hit = False

  # Initialize first individuals 
  shape = theta
  scale = nubar / shape
  nu0 = np.random.gamma(shape, scale, ninit)
  nu[:ninit] = nu0
  on[:ninit] = np.random.poisson(nu0)
  gen[:ninit] = 0

  # Simulation loop
  i = 0  # index of current focal idvl
  j = ninit  # index of next empty position in the arrays

  while True:
    
    # Stop if no more created individuals to process 
    # (or safety cap)
    if i >= j or i >= boundary:
      break
    
    # Get data for focal idvl
    nu_focal = nu[i]  # irn focal
    on_focal = int(on[i])  # offspring number focal

    # Write these offspring
    if on_focal > 0:

      remaining = boundary - j
      if remaining <= 0:
        boundary_hit = True
        break

      truncated = False
      if on_focal > remaining:
        on_focal = remaining
        on[i] = on_focal
        truncated = True

      # generate the offspring that fit
      nup_focal = nu_focal
      mu = (1 - delta) * nubar + delta * nup_focal
      shape = theta
      scale = mu / shape

      nu_secondary = np.random.gamma(shape, scale, on_focal)
      on_secondary = np.random.poisson(nu_secondary)

      x = np.arange(j, j + on_focal)

      parent[x] = i
      nup[x] = nup_focal
      nu[x] = nu_secondary
      onp[x] = on_focal
      on[x] = on_secondary
      gen[x] = gen[i] + 1

      j += on_focal

      if truncated:
        boundary_hit = True
        break
      
    # Advance to next focal individual
    i += 1
        
  # Calculate average individual reproductive number per generation
  # and also average parent individual reproductive number per generation
  unique_gen = np.unique(gen[~np.isnan(gen)])
  r   = np.array([np.mean(nu[gen == g]) for g in unique_gen])
  r_p = np.array([np.mean(nup[gen == g]) for g in unique_gen])
  r   = r[:ngen] if len(r) > ngen else np.pad(r, (0, ngen - len(r)), constant_values = np.nan)
  r_p = r_p[:ngen] if len(r_p) > ngen else np.pad(r_p, (0, ngen - len(r_p)), constant_values = np.nan)    
      
      
      
  # Calculate the size of each generation
  unique_gen = np.unique(gen[~np.isnan(gen)])
  gen_size = np.array([np.sum(gen == g) for g in unique_gen], dtype=float)
  gen_size = gen_size[:ngen] if len(gen_size) > ngen else np.pad(gen_size, (0, ngen - len(gen_size)), constant_values = np.nan)


  # Assemble data
  data = pd.DataFrame({
      'parent': parent,
      'nup': nup,
      'nu': nu,
      'onp': onp,
      'on': on,
      'gen': gen
  })
  data = data.iloc[:j].reset_index(drop=True)


  # Return
  final_size = j
  final_gen = int(np.nanmax(gen)) 

  out = {
      'data': data,
      'final_size': final_size,
      'boundary_hit': boundary_hit,
      'r': r,
      'r_p': r_p,
      'gen_size': gen_size,
      'final_gen': final_gen
  }
  return out