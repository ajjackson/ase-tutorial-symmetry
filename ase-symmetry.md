---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.1'
      jupytext_version: 1.2.4
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Symmetry in ASE


## Using symmetry to build structures


## Bravais lattices
The simplest crystal structures are defined entirely by their periodicity; the position of the atom within the unit cell is unimportant.

Here we consider metallic silver: the list of available Bravais lattices in ASE is available [here](https://wiki.fysik.dtu.dk/ase/ase/lattice.html#available-crystal-lattices).

```python
import ase.lattice

ag_lattice = ase.lattice.FCC(a=4.09) # Face-centered cubic
print("Made a lattice for Ag: ", type(ag_lattice))

assert isinstance(ag_lattice, ase.lattice.BravaisLattice)
```

This is an instance of the FCC class which is subclassed from BravaisLattice. (We don't work with BravaisLattice directly.) We can use this instance to construct the unit cell for our Atoms object.

```python
from ase import Atoms
from ase.visualize import view
ag_cell = ag_lattice.tocell()
ag_atoms = Atoms('Ag', cell=ag_cell, pbc=True)
view(ag_atoms)
print(ag_atoms.cell)
```

Note that this has constructed a *primitive* cell of Ag rather than the cubic cell with 4.09Å sides.


> It's not obvious what the lattice generating functions will be named.
> We could look in the documentation or the source code, but from an iPython terminal
> or Jupyter notebook there is another convenient option: tab-completion.
> Type `ase.lattice.` and hit the tab key a few times to open a list of available options.
> These are all the classes, functions and properties available in the `ase.lattice` module.

> *Find the tetragonal lattice and build a tetragonal unit cell. Have a look in the viewer to see if it makes sense.*


## Space groups


In crystallography structures are typically defined by a combination of
- Bravais lattice
- Symmetry operations
- Atomic basis

All the possible combinations of Bravais lattice types and symmetry operations have been determined by group theory and are assigned to crystallogaphic space groups. There are a few notation schemes, and arguably the easiest to work with are the sequential numbers 1-230 given in the *International Tables for Crystallography, vol. A*.


Examining a public-domain .cif file from Crystallography Open Database:

```bash
cat quartz.cif
```

<!-- #region -->
Working from the bottom to the top of the file, we see in the last two lines that coordinates are only provided for two atoms in the cell. The "loop_" lines above define that format for these coordinates -- CIF is not the most concise format for structure data! Above these lines we see a loop over symmetry operations:
```
loop_
_symmetry_equiv_pos_as_xyz
x,y,z
-y,x-y,z+1/3
-x+y,-x,z+2/3
y,x,-z
-x,-x+y,-z+1/3
x-y,-y,-z+2/3
```

`x,y,z` is the identity operation and `y,x,-z` is a rotation/reflection operation. The other operations include translations in the z axis - given that we only have one Si and one O in the basis, these *must* create some extra atoms.

At the top of the file are lattice parameters and metadata, including several definitions of the spacegroup.
```
_space_group_IT_number           152
_symmetry_space_group_name_Hall  'P 31 2"'
_symmetry_space_group_name_H-M   'P 31 2 1'
```
<!-- #endregion -->

```python
import ase.spacegroup
quartz = ase.spacegroup.crystal(symbols=['O', 'Si'],
                                basis=[[0.413, 0.2711, 0.2172],
                                       [0.4673, 0, 0.3333]],
                                spacegroup=152,
                                cellpar=[4.9019, 4.9019, 5.3988, 90, 90, 120])
view(quartz)
```

Hopefully you can convince yourself in the viewer that this is a network of tetrahedral units, consisting of Si surrounded by O. It may be helpful in the ASE viewer to use the "repeats" option in the View menu to expand the viewing region. Alternatively, we can dump out a structure with e.g.

```python
quartz.write('POSCAR', vasp5=True)
```

for viewing in your favourite atomistic structure viewer. (Mine is VESTA.) You can also open up the original .cif file and see if the results agree.

You don't *have* to make structures from .cif files this way as ASE has CIF support built-in; we could have saved ourselves a lot of trouble with

```python
import ase.io
atoms = ase.io.read('quartz.cif')
```

but it's not unusual to find structures specified this way in research papers or textbooks.


### Getting the symmetry operations
How did ASE do that?

ASE has a data structure for spacegroups that includes information about the symmetry operations, loaded from a data table.

```python
spg = ase.spacegroup.Spacegroup(152)
print("Is a-quartz centrosymmetric? {}".format("yes" if spg.centrosymmetric else "no"))
print("Is the a-quartz lattice primitive? {}".format("yes" if spg.lattice == 'P' else "no"))
print("What symbol is this spacegroup? {}".format(spg.symbol))
```

Perhaps most importantly the actual symmetry operations are made available a list or iterator.

```python
for i, (rot, trans) in enumerate(spg.get_symop()):
    print("Symmetry operation #{}".format(i + 1))
    print("Rotation matrix: ", rot[0,:])
    print("                 ", rot[1,:])
    print("                 ", rot[2,:])
    print("Translation vector: ", trans)
    print()
```

Hopefully these operations look familiar from the CIF file!

> *Have a look at the symmetry operations from some other spacegroups. Can you identify rotation, reflection and glide operations?*


To generate a structure with these operations, ASE uses the symops in a method `Spacegroup.equivalent_sites`. This takes positions in fractional coordinates and applies the symmetry operations to derive a set of images. If we provide just the Si atom from the .cif file, we get our original site back, with two other images - this is consistent with the cell of three formula units.

```python
spg.equivalent_sites([0.4673, 0, 0.3333])
```

> *What does the [0, 0, 0] mean? This is not part of the equivalent_sites array. Have a look at the docstring...*
>
> To examine detailed help from a regular python terminal you could use `help(ase.spacegroup.Spacegroup)`.
> To inspect the manually-written docstring part directly you could print `ase.spacegroup.Spacegroup.__doc__`.
> From an iPython terminal or Jupyter notebook the most convenient option is the special syntax `ase.spacegroup.Spacegroup?`.
> Most functions in ASE should have docstrings, so adding a `?` to a half-written line of code can be very helpful in the heat of the moment.
