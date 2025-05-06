#!/usr/bin/env python3

import numpy as np

tipsy_header_type = np.dtype([('time', '>f8'),('N', '>u4'), ('Dims', '>u4'), ('Ngas', '>u4'), ('Ndark', '>u4'), ('Nstar', '>u4'), ('pad', '>u4')])
gas_type  = np.dtype([('mass','>f4'),('x','>f4'),('y','>f4'),('z','>f4'),('vx','>f4'),('vy','>f4'),('vz','>f4'),('rho','>f4'),('temp','>f4'),('hsmooth','>f4'),('metals','>f4'),('phi','>f4')])

def main():
    dKpcUnit = 2.06701e-13
    dMsolUnit = 4.80438e-08
    N = 100e6
    R = 1
    V = 4.0 * np.pi / 3.0 * R**3.0
    a_cell = (4.0 * np.pi / (N/2.0 * 3.0) * R**3.0)**(1.0/3.0)

    m = 62.5476 / N
   
    a = -R*1.1
    b = R*1.1

    # Create initial grid
    Nx = int(np.ceil((b-a)/a_cell))
    x = np.linspace(-0.5*Nx*a_cell, 0.5*Nx*a_cell, Nx+1)
    x_1, y_1, z_1 = np.meshgrid(x, x, x, indexing='ij')
    # Create offset grid
    x_2 = x_1 + 0.5 * a_cell
    y_2 = y_1 + 0.5 * a_cell
    z_2 = z_1 + 0.5 * a_cell
    # Remove particles outside of sphere
    r_1 = np.sqrt(x_1**2.0 + y_1**2.0 + z_1**2.0)
    outside = r_1 >= R
    x_1 = x_1[~outside]
    y_1 = y_1[~outside]
    z_1 = z_1[~outside]
    r_2 = np.sqrt(x_2**2.0 + y_2**2.0 + z_2**2.0)
    outside = r_2 >= R
    x_2 = x_2[~outside]
    y_2 = y_2[~outside]
    z_2 = z_2[~outside]

    # Create gas structure for initial grid
    gas_1 = np.zeros(np.shape(x_1)[0],dtype=gas_type)
    gas_1['mass'] = m
    gas_1['x'] = x_1
    gas_1['y'] = y_1
    gas_1['z'] = z_1
    gas_1['vx'] = 0.0
    gas_1['vy'] = 0.0
    gas_1['vz'] = 0.0
    gas_1['rho'] = 0.0
    gas_1['temp'] = 0.0
    gas_1['hsmooth'] = 0.0
    gas_1['metals'] = 101.0
    gas_1['phi'] = 0.0
    # Create gas structure for offset grid
    gas_2 = np.zeros(np.shape(x_2)[0],dtype=gas_type)
    gas_2['mass'] = m
    gas_2['x'] = x_2
    gas_2['y'] = y_2
    gas_2['z'] = z_2
    gas_2['vx'] = 0.0
    gas_2['vy'] = 0.0
    gas_2['vz'] = 0.0
    gas_2['rho'] = 0.0
    gas_2['temp'] = 0.0
    gas_2['hsmooth'] = 0.0
    gas_2['metals'] = 102.0
    gas_2['phi'] = 0.0

    # Generate IC for fully segregated
    gas = np.append(gas_1,gas_2)
    r = np.sqrt(gas['x']*gas['x'] + gas['y']*gas['y'] + gas['z']*gas['z'])
    R_core = (1.0 / 3.0)**(1.0/3.0)
    gas['metals'][r<R_core] = 1.0
    gas['metals'][r>=R_core] = 0.0
    # Write tipsy file
    header = np.zeros(1,dtype=tipsy_header_type)
    header['time'] = 0.0
    header['N'] = np.shape(gas)[0]
    header['Dims'] = 3
    header['Ngas'] = np.shape(gas)[0]
    header['Ndark'] = 0
    header['Nstar'] = 0
    header['pad'] = 0
    newfile = open('01_fully_segregated.std','wb')
    newfile.write(np.array(header,dtype=tipsy_header_type).tobytes())
    newfile.write(np.array(gas,dtype=gas_type).tobytes())
    newfile.close()
    
    # Generate IC for central sphere even/odd
    gas = np.append(gas_1,gas_2)
    r = np.sqrt(gas['x']*gas['x'] + gas['y']*gas['y'] + gas['z']*gas['z'])
    R_core = (2.0 / 3.0)**(1.0/3.0)
    gas['metals'][np.logical_and(r<R_core,gas['metals']==102.0)] = 1.0
    gas['metals'][np.logical_and(r<R_core,gas['metals']==101.0)] = 0.0
    gas['metals'][r>=R_core] = 0.0
    # Write tipsy file
    header = np.zeros(1,dtype=tipsy_header_type)
    header['time'] = 0.0
    header['N'] = np.shape(gas)[0]
    header['Dims'] = 3
    header['Ngas'] = np.shape(gas)[0]
    header['Ndark'] = 0
    header['Nstar'] = 0
    header['pad'] = 0
    newfile = open('02_central_sphere_even_odd.std','wb')
    newfile.write(np.array(header,dtype=tipsy_header_type).tobytes())
    newfile.write(np.array(gas,dtype=gas_type).tobytes())
    newfile.close()
    
    # Generate IC for central sphere random
    gas = np.append(gas_1,gas_2)
    r = np.sqrt(gas['x']*gas['x'] + gas['y']*gas['y'] + gas['z']*gas['z'])
    R_core = (2.0 / 3.0)**(1.0/3.0)
    gas['metals'][r<R_core] = np.random.randint(0, high=2, size=np.count_nonzero(gas['metals'][r<R_core])) < 0.5
    gas['metals'][r>=R_core] = 0.0
    # Write tipsy file
    header = np.zeros(1,dtype=tipsy_header_type)
    header['time'] = 0.0
    header['N'] = np.shape(gas)[0]
    header['Dims'] = 3
    header['Ngas'] = np.shape(gas)[0]
    header['Ndark'] = 0
    header['Nstar'] = 0
    header['pad'] = 0
    newfile = open('03_central_sphere_random.std','wb')
    newfile.write(np.array(header,dtype=tipsy_header_type).tobytes())
    newfile.write(np.array(gas,dtype=gas_type).tobytes())
    newfile.close()
    
    # Generate IC for full sphere random
    gas = np.append(gas_1,gas_2)
    r = np.sqrt(gas['x']*gas['x'] + gas['y']*gas['y'] + gas['z']*gas['z'])
    R_core = (2.0 / 3.0)**(1.0/3.0)
    gas['metals'] = np.random.randint(0, high=3, size=np.count_nonzero(gas['metals'])) < 1.0/3.0
    # Write tipsy file
    header = np.zeros(1,dtype=tipsy_header_type)
    header['time'] = 0.0
    header['N'] = np.shape(gas)[0]
    header['Dims'] = 3
    header['Ngas'] = np.shape(gas)[0]
    header['Ndark'] = 0
    header['Nstar'] = 0
    header['pad'] = 0
    newfile = open('04_full_sphere_random.std','wb')
    newfile.write(np.array(header,dtype=tipsy_header_type).tobytes())
    newfile.write(np.array(gas,dtype=gas_type).tobytes())
    newfile.close()
    
    # Create initial grid
    Nx = int(np.ceil((b-a)/a_cell))
    x = np.linspace(-0.5*Nx*a_cell, 0.5*Nx*a_cell, Nx+1)
    x_1, y_1, z_1 = np.meshgrid(x, x, x, indexing='ij')
    # Create offset grid
    x_2 = x_1 + 0.5 * a_cell
    y_2 = y_1 + 0.5 * a_cell
    z_2 = z_1 + 0.5 * a_cell
    spread = 0.5
    x_1 += 2.0 * spread * a_cell * np.random.random_sample(np.shape(x_1)) - spread * a_cell
    y_1 += 2.0 * spread * a_cell * np.random.random_sample(np.shape(y_1)) - spread * a_cell
    z_1 += 2.0 * spread * a_cell * np.random.random_sample(np.shape(z_1)) - spread * a_cell
    x_2 += 2.0 * spread * a_cell * np.random.random_sample(np.shape(x_2)) - spread * a_cell
    y_2 += 2.0 * spread * a_cell * np.random.random_sample(np.shape(y_2)) - spread * a_cell
    z_2 += 2.0 * spread * a_cell * np.random.random_sample(np.shape(z_2)) - spread * a_cell
    # Remove particles outside of sphere
    r_1 = np.sqrt(x_1**2.0 + y_1**2.0 + z_1**2.0)
    outside = r_1 >= R
    x_1 = x_1[~outside]
    y_1 = y_1[~outside]
    z_1 = z_1[~outside]
    r_2 = np.sqrt(x_2**2.0 + y_2**2.0 + z_2**2.0)
    outside = r_2 >= R
    x_2 = x_2[~outside]
    y_2 = y_2[~outside]
    z_2 = z_2[~outside]
    
    # Create gas structure for initial grid
    gas_1 = np.zeros(np.shape(x_1)[0],dtype=gas_type)
    gas_1['mass'] = m
    gas_1['x'] = x_1
    gas_1['y'] = y_1
    gas_1['z'] = z_1
    gas_1['vx'] = 0.0
    gas_1['vy'] = 0.0
    gas_1['vz'] = 0.0
    gas_1['rho'] = 0.0
    gas_1['temp'] = 0.0
    gas_1['hsmooth'] = 0.0
    gas_1['metals'] = 101.0
    gas_1['phi'] = 0.0
    # Create gas structure for offset grid
    gas_2 = np.zeros(np.shape(x_2)[0],dtype=gas_type)
    gas_2['mass'] = m
    gas_2['x'] = x_2
    gas_2['y'] = y_2
    gas_2['z'] = z_2
    gas_2['vx'] = 0.0
    gas_2['vy'] = 0.0
    gas_2['vz'] = 0.0
    gas_2['rho'] = 0.0
    gas_2['temp'] = 0.0
    gas_2['hsmooth'] = 0.0
    gas_2['metals'] = 102.0
    gas_2['phi'] = 0.0

    # Generate IC for fully segregated
    gas = np.append(gas_1,gas_2)
    r = np.sqrt(gas['x']*gas['x'] + gas['y']*gas['y'] + gas['z']*gas['z'])
    R_core = (1.0 / 3.0)**(1.0/3.0)
    gas['metals'][r<R_core] = 1.0
    gas['metals'][r>=R_core] = 0.0
    # Write tipsy file
    header = np.zeros(1,dtype=tipsy_header_type)
    header['time'] = 0.0
    header['N'] = np.shape(gas)[0]
    header['Dims'] = 3
    header['Ngas'] = np.shape(gas)[0]
    header['Ndark'] = 0
    header['Nstar'] = 0
    header['pad'] = 0
    newfile = open('05_fully_segregated_moved.std','wb')
    newfile.write(np.array(header,dtype=tipsy_header_type).tobytes())
    newfile.write(np.array(gas,dtype=gas_type).tobytes())
    newfile.close()
    
    # Generate IC for central sphere even/odd
    gas = np.append(gas_1,gas_2)
    r = np.sqrt(gas['x']*gas['x'] + gas['y']*gas['y'] + gas['z']*gas['z'])
    R_core = (2.0 / 3.0)**(1.0/3.0)
    gas['metals'][np.logical_and(r<R_core,gas['metals']==102.0)] = 1.0
    gas['metals'][np.logical_and(r<R_core,gas['metals']==101.0)] = 0.0
    gas['metals'][r>=R_core] = 0.0
    # Write tipsy file
    header = np.zeros(1,dtype=tipsy_header_type)
    header['time'] = 0.0
    header['N'] = np.shape(gas)[0]
    header['Dims'] = 3
    header['Ngas'] = np.shape(gas)[0]
    header['Ndark'] = 0
    header['Nstar'] = 0
    header['pad'] = 0
    newfile = open('06_central_sphere_even_odd_moved.std','wb')
    newfile.write(np.array(header,dtype=tipsy_header_type).tobytes())
    newfile.write(np.array(gas,dtype=gas_type).tobytes())
    newfile.close()
    
    # Generate IC for central sphere random
    gas = np.append(gas_1,gas_2)
    r = np.sqrt(gas['x']*gas['x'] + gas['y']*gas['y'] + gas['z']*gas['z'])
    R_core = (2.0 / 3.0)**(1.0/3.0)
    gas['metals'][r<R_core] = np.random.randint(0, high=2, size=np.count_nonzero(gas['metals'][r<R_core])) < 0.5
    gas['metals'][r>=R_core] = 0.0
    # Write tipsy file
    header = np.zeros(1,dtype=tipsy_header_type)
    header['time'] = 0.0
    header['N'] = np.shape(gas)[0]
    header['Dims'] = 3
    header['Ngas'] = np.shape(gas)[0]
    header['Ndark'] = 0
    header['Nstar'] = 0
    header['pad'] = 0
    newfile = open('07_central_sphere_random_moved.std','wb')
    newfile.write(np.array(header,dtype=tipsy_header_type).tobytes())
    newfile.write(np.array(gas,dtype=gas_type).tobytes())
    newfile.close()
    
    # Generate IC for full sphere random
    gas = np.append(gas_1,gas_2)
    r = np.sqrt(gas['x']*gas['x'] + gas['y']*gas['y'] + gas['z']*gas['z'])
    R_core = (2.0 / 3.0)**(1.0/3.0)
    gas['metals'] = np.random.randint(0, high=3, size=np.count_nonzero(gas['metals'])) < 1.0/3.0
    # Write tipsy file
    header = np.zeros(1,dtype=tipsy_header_type)
    header['time'] = 0.0
    header['N'] = np.shape(gas)[0]
    header['Dims'] = 3
    header['Ngas'] = np.shape(gas)[0]
    header['Ndark'] = 0
    header['Nstar'] = 0
    header['pad'] = 0
    newfile = open('08_full_sphere_random_moved.std','wb')
    newfile.write(np.array(header,dtype=tipsy_header_type).tobytes())
    newfile.write(np.array(gas,dtype=gas_type).tobytes())
    newfile.close()
    

if __name__ == '__main__':
    main()
