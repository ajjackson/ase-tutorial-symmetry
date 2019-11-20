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

# Crystal symmetry in ASE

This notebook provides an introduction to some of the crystal symmetry features in ASE and spglib.

It is intended as a self-guided tour, with small problems set in **bold text**. A few example systems are used, and readers are strongly encouraged to try out these tools with systems they are interested in, especially those with non-primitive unit cells and where data is available from geometry optimisation or structure prediction.

> Notes written in quote blocks give tips that are related more to Python, IPython and Jupyter features.
> Some familiarity with Python is assumed, but these tips should be helpful for completing the problems.

The tutorial will first look at how we can apply known symmetry information when creating a system, then move on to analysing the symmetry of an existing structure.


## Bravais lattices
The simplest crystal structures are defined entirely by their periodicity; the position of the atom within the unit cell is unimportant.

Here we consider metallic silver: the list of available Bravais lattices in ASE is available [here](https://wiki.fysik.dtu.dk/ase/ase/lattice.html#available-crystal-lattices).

```python
import ase.lattice

ag_lattice = ase.lattice.FCC(a=4.09) # Face-centered cubic
print("Made a lattice for Ag: ", type(ag_lattice))

assert isinstance(ag_lattice, ase.lattice.BravaisLattice)
```

This is an instance of the FCC class which is derived from a parent class *BravaisLattice*. We don't work with BravaisLattice directly, but it provides features to its "children". Here we have set up an instance of FCC with lattice parameters matching metallic silver. To make an Atoms object, we convert to a unit cell and then to Atoms.
> In Python, we can compare types with `type(x) == type(y)` but this will only succeed if they are exactly the same.
> It is often more robust to use `isinstance()`, which will return True if the object type matches _or inherits from_ the compared type.

```python
from ase import Atoms
ag_cell = ag_lattice.tocell()
ag_atoms = Atoms('Ag', cell=ag_cell, pbc=True)
print(ag_atoms.cell)
```

Note that _Atoms_ needs to know what element we are using, while the _Cell_ does not. Because there is one atom and it doesn't matter where it goes, we didn't have to specify any positions.

This operation has constructed a *primitive* cell of Ag rather than the cubic cell with 4.09Å sides.
We can view the structure with ASE's inbuilt viewer. By default the `ase.visualize.view` command will create a pop-up window with the structure in it (you may need to move your browser window to spot it!)

```python
from ase.visualize import view
view(ag_atoms)
```

A more notebook-friendly option is to plot with matplotlib. (Interactive notebook viewers are also avilable; see the [ASE docs](https://wiki.fysik.dtu.dk/ase/ase/visualize/visualize.html#viewer-for-jupyter-notebooks).)

```python
%matplotlib notebook
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
plot_atoms(ag_atoms)
```

**Now try to build a tetragonal unit cell. Visualise it to see if it makes sense.**

> It's not obvious what the lattice generating functions will be named.
> We could look in the documentation or the source code, but from an iPython terminal
> or Jupyter notebook there is another convenient option: tab-completion.
> Type `ase.lattice.` and hit the tab key a few times to open a list of available options.
> These are all the classes, functions and properties available in the `ase.lattice` module.



## Space groups


In crystallography, material structures are typically defined by a combination of
- Bravais lattice
- Symmetry operations
- Atomic basis

All the possible combinations of Bravais lattice types and symmetry operations have been determined by group theory and are assigned to crystallogaphic space groups. There are a few notation schemes, and arguably the easiest to work with are the sequential numbers 1-230 given in the *International Tables for Crystallography, vol. A*.


Let's examine a public-domain .cif file from Crystallography Open Database. Put this file, downloaded from www.crystallography.net/cod/cif/1/52/68/1526860.cif, in the current working directory a rename it to *quartz.cif*.

> From a Linux or Mac system it's possible to make a call to the Bash shell to display the file, even though this is primarily a Python notebook, with the "cell magic" `%%bash`. You can also view this file with a text editor.

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

In the `ase.spacegroup` module we find the command `ase.spacegroup.crystal` which allows us to build an _Atoms_ object using this crystallographic information.

> This may not be the most obvious place as a user to look for this function: other high-symmetry "structure builder" functions live in `ase.build`. However, this is the module that contains the underlying symmetry tools, as we'll see in the next part of the tutorial.

```python
import ase.spacegroup
quartz = ase.spacegroup.crystal(symbols=['O', 'Si'],
                                basis=[[0.413, 0.2711, 0.2172],
                                       [0.4673, 0, 0.3333]],
                                spacegroup=152,
                                cellpar=[4.9019, 4.9019, 5.3988, 90, 90, 120])
fig, ax = plt.subplots()
plot_atoms(quartz)
```

**Convince yourself that this structure is sensible**

> Hopefully you can find out that this is a network of tetrahedral units, consisting of Si surrounded by O. It may be helpful to open the 3D viewer (using `view()`) and use the "repeats" option in the View menu to expand the viewing region. Alternatively, we can dump out a structure with e.g. `quartz.write('POSCAR', vasp5=True)` for viewing in your favourite atomistic structure viewer (e.g. VESTA).


You can also open up the original .cif file and see if the results agree.

You don't *have* to make structures from .cif files this way. ASE has CIF support built-in, and we could have saved ourselves a lot of trouble with

```python
import ase.io
atoms = ase.io.read('quartz.cif')
```

but it's not unusual to find structures specified this way in research papers or textbooks.


**Examine a .cif file for a high-symmetry material of interest to your research. How many symmetry operations are there? What happens when you construct the cell in ASE?**


### Getting the symmetry operations
How did ASE do that?

ASE has a data structure for spacegroups that includes information about the symmetry operations, loaded from a data table.

```python
spg = ase.spacegroup.Spacegroup(152)
print("Is a-quartz centrosymmetric? {}".format("yes" if spg.centrosymmetric else "no"))
print("Is the a-quartz lattice primitive? {}".format("yes" if spg.lattice == 'P' else "no"))
print("What symbol is this spacegroup? {}".format(spg.symbol))
```

This Class provides access to the actual symmetry operations as a list or iterator.

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

**Have a look at the symmetry operations from some other spacegroups. Can you identify rotation, reflection and glide operations?**

To generate a structure with these operations, ASE uses the symops in a method `Spacegroup.equivalent_sites`. This takes positions in fractional coordinates and applies the symmetry operations to derive a set of images. If we provide just the Si atom from the .cif file, we get our original site back, with two other images - this is consistent with the cell of three formula units.

```python
spg.equivalent_sites([0.4673, 0, 0.3333])
```

**What does the [0, 0, 0] mean? This is not part of the equivalent_sites array. Have a look at the docstring...**
>
> To examine detailed help from a regular python terminal you could use `help(ase.spacegroup.Spacegroup)`.
> To inspect the manually-written docstring part directly you could print `ase.spacegroup.Spacegroup.__doc__`.
> From an iPython terminal or Jupyter notebook the most convenient option is the special syntax `ase.spacegroup.Spacegroup?`.
> Most functions in ASE should have docstrings, so adding a `?` to a half-written line of code can be very helpful in the heat of the moment.

**What happens if you include redundant sites in the input to *equivalent_sites()*?**


## Analysing symmetry in existing structures

Many atomistic codes are "smart" about how they use symmetry. They can identify symmetry and use this to:

- reduce the amount of information to be calculated;
- reduce the size of data files;
- impose constraints on calculated structures and fields.

However, there are several reasons reasearchers _should_ be aware (and possibly worried) about the symmetry of their system, including:

- working with primitive/conventional unit cells as appropriate;
- identifying unique sites for substitutions and defects;
- identifying when symmetry was broken by numerical errors;
- avoiding the imposition of excessive symmetry on top of the _magnetic_ spacegroup.


### Identifying the spacegroup

Going back to the quartz cell that was imported from a .cif file, we can have a look at some metadata stored in the `info` attribute of *Atoms*.

```python
atoms.info
```

The spacegroup was actually set as a tag when it was imported. Let's break the symmetry slightly by making some random changes to the atom positions.

```python
atoms.rattle()
atoms.positions
```

This didn't change the tag, so to get the correct symmetry we need to analyse the structure.

```python
ase.spacegroup.get_spacegroup(atoms)
```

We now have spacegroup 1, with zero symmetry. That's an accurate assessment of the system, but there are plenty of cases where we have a structure that isn't quite _perfect_ and we'd like to know what high-symmetry structure is nearby. To do this we can increase the tolerance used to determine if atoms are "in the same position".

```python
ase.spacegroup.get_spacegroup(atoms, symprec=1e-2)
```

This useful function is actually provided by [spglib](https://atztogo.github.io/spglib), an external library written in C. The performance benefits of C are quite helpful here, as determining the spacegroup involves applying many symmetry operations and analysing the results. Spglib includes a Python API and the integration into ASE is not too complicated. We can also access the same function directly

```python
from spglib import get_spacegroup
get_spacegroup((atoms.cell, atoms.get_scaled_positions(), atoms.numbers),
               symprec=1e-2)
```

It's a bit more cumbersome compared to calling directly from ASE because spglib doesn't know how to read an Atoms object, so we extract the necessary parts and put them in a `(tuple)`.
**When might it be useful to call this function directly? Look at the docstring...**


### Cleaning up a structure with standardize_cell

There are a few more functions in the *spglib interface* that are not directly provided in ASE.

Suppose your collaborator has asked you to calculate some properties and has sent you the structure file *mystery.cell* (included in the directory with this notebook).

```python
mystery_atoms = ase.io.read('mystery.cell')
ase.spacegroup.get_spacegroup(mystery_atoms)
```

Well, you _could_ assume the structure is already perfect and we are simply dealing with a low-symmetry system.

**Visualise the structure. Does this look low-symmetry to you?**

It's not always easy to judge by eye. We don't know what spacegroup the material is _supposed_ to have. A useful procedure in this scenario is to look at how sensitive the spacegroup is to the distance threshold.

```python
for symprec in (1e-5, 5e-4, 1e-4, 1e-3, 5e-3, 1e-2, 1e-1):
    print("Threshold: {}  Spacegroup: {}".format(symprec,
                                                 ase.spacegroup.get_spacegroup(mystery_atoms,
                                                                               symprec=symprec).symbol))
```

So there is an I-4 system within a small atomic movement. It's possible that the symmetry-breaking is real, of course, but it would be nice to do some calculations with the high-symmetry system and get those efficiency improvements.

```python
from spglib import standardize_cell
```

`standardize_cell` replaces the older `find_primitive` and `refine_cell` functions, and is generally recommended for this purpose.

```python
mystery_cell = (mystery_atoms.cell, mystery_atoms.get_scaled_positions(), mystery_atoms.numbers)
new_cell = standardize_cell(mystery_cell, to_primitive=False, symprec=5e-3)
new_cell
```

Note that an spglib _cell_ is more analogous to _Atoms_ in ASE than to a _Cell_ in ASE. To convert this to an Atoms object we can feed the data into an Atoms constructor. This is made especially clear using Python's variable-unpacking syntax `x, y, z = a`.

```python
new_unit_cell, new_scaled_positions, new_numbers = new_cell
new_atoms = Atoms(new_numbers, cell=new_unit_cell, scaled_positions=new_scaled_positions)
fig, ax = plt.subplots()
plot_atoms(new_atoms, rotation='90x,90y')
```

**Write this structure to a .cell file and compare it with the input mystery.cell file**

This cell contains more atoms than the original data file. For many purposes this _conventional_ unit cell is preferred, but generally it will be more efficient to perform further calculations using the smaller _primitive_ cell.

**Generate a primitive cell instead and write this to your preferred geometry file format.**
> `spglib.standardize_cell` contains what you need to do this.


**Final exercise:**
The "mystery" system is a photovoltaic absorber material which has been found to consistently under-perform compared to its theoretical efficiency limit. This might be due to its defect chemistry. The formation energies of isolated vacancies can be computed by constructing a supercell and removing appropriate atoms. How many inequivalent atom sites do we need to consider in this material?

> Have a look at the available methods to Spacegroup and see which ones might be useful. There are several options!

(This kind of study has been performed a few times, and recent work with hybrid DFT suggests that sulfur vacancies do play a significant role in this system.)
