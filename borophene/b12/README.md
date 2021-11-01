## Initial Structures

b12_11.in (1x1 structure for conventional cell of b12)

for the b12 phase  vc_realxed structure (only cell parameters relaxed and atomic positions are in their high symmetry positions ) 

cell parameters are

```bash
celldm(1)=9.42608683764, celldm(2)=0.5922, celldm(3)=1.2d0,
```

in other format 

```bash
4.98806835 0.0000 0.0
0.00000000 2.9539 0.0
0.00000000 0.0000 6.0
```



Then we made a 2x 2 supercell and optimized its lattice constants in xy direction

```bash
&cell
cell_dofree='xy'
```

and finally we find bellow parameters for 2x2 supercell host 

```bash
celldm(1)=19.16176391819488, celldm(2)=0.57671, celldm(3)=0.700000000
```



OR  in angstrom

```bash
10.139964 0.000000 0.00
0.0000000 5.847908 0.00
0.0000000 0.000000 7.09
```





