import numpy as np

### from SUEWS

top_flux_dn_sw = 200
top_flux_dn_direct_sw = 300
top_flux_dn_lw = 300

#### bc_out:
    
sw_albedo = 0.11325798862494915
sw_albedo_dir = 0.13138550212367153
lw_emissivity = 0.94721882875278507
lw_emission = 160.38901583274799

#### lw_flux (in-out of facet or down-up): 
    
ground_dn = 29.522044579855191

top_dn = 300

roof_in = np.array([0,0,150])
wall_in = np.array([107.58466809856942,102.02507757929145,109.16598121932485])

clear_air_abs = np.array([0,0,0])
veg_abs = np.array([126.5167503042217,70.425781939607603,137.22542562810378])
veg_air_abs = np.array([0,0,0])

wall_net = np.array([-34.945566668442304,-39.949198003241527,-33.522384897463866])
ground_net = -105.20192597412712
roof_net = np.array([0,0,3.2282310315848974])

top_net = 123.77663279308751

### lw top_net
print("\nLW")

sum_clear_air_abs = np.sum(clear_air_abs)
print("\nsum_clear_air_abs:",sum_clear_air_abs)
sum_veg_abs = np.sum(veg_abs)
print("sum_veg_abs:",sum_veg_abs)
sum_veg_air_abs = np.sum(veg_air_abs)
print("sum_veg_air_abs:",sum_veg_air_abs)
sum_wall_net = np.sum(wall_net)
print("sum_wall_net:",sum_wall_net)
sum_ground_net = np.sum(ground_net)
print("sum_ground_net:",sum_ground_net)
sum_roof_net = np.sum(roof_net)
print("sum_roof_net:",sum_roof_net)

top_net = sum_roof_net + sum_ground_net + sum_wall_net \
            + sum_clear_air_abs + sum_veg_abs + sum_veg_air_abs
print("\nlw top_net from facets, air and veg:",top_net)

top_net = top_flux_dn_lw - (top_flux_dn_lw*(1-lw_emissivity) + lw_emission)  #  top_flux_dn_lw*lw_emissivity - lw_emission
print("lw top_net from bc_out emissivity and emission:",top_net)

#### sw_flux (in-out of facet or down-up):
    
ground_dn = -0.44753600582403369
ground_dn_dir = 0

top_dn = 200
top_dn_dir = 300

roof_in = np.array([0,0,100])
wall_in = np.array([-0.43954298655747076,-1.6098738909684078,27.574735360648152])

clear_air_abs = np.array([-0.022084275784198578,-0.080706822785717655,0.74854740883402271])
veg_abs = np.array([-1.0472198505274308,-4.5171661162312313,73.666721895426718])
veg_air_abs = np.array([-0.047666557397738681,-0.20560893478672299,3.3531058695679601])

wall_net = np.array([-0.35163438793603646,-1.2878991079769284,22.059788206339409])
ground_net = -0.35802880332546577
roof_net = np.array([0,0,79.999999701976776])

top_net = 171.91014822539347

### sw top_net
print("\nSW")

sum_clear_air_abs = np.sum(clear_air_abs)
print("\nsum_clear_air_abs:",sum_clear_air_abs)
sum_veg_abs = np.sum(veg_abs)
print("sum_veg_abs:",sum_veg_abs)
sum_veg_air_abs = np.sum(veg_air_abs)
print("sum_veg_air_abs:",sum_veg_air_abs)
sum_wall_net = np.sum(wall_net)
print("sum_wall_net:",sum_wall_net)
sum_ground_net = np.sum(ground_net)
print("sum_ground_net:",sum_ground_net)
sum_roof_net = np.sum(roof_net)
print("sum_roof_net:",sum_roof_net)

top_net = sum_roof_net + sum_ground_net + sum_wall_net \
            + sum_clear_air_abs + sum_veg_abs + sum_veg_air_abs
print("\nsw top_net from facets, air and veg:",top_net)

top_dn_diff = top_dn-top_dn_dir

top_net = top_flux_dn_direct_sw*(1-sw_albedo_dir) + top_dn_diff*(1-sw_albedo)
print("sw top_net from bc_out albedos:",top_net)