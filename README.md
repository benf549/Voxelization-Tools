# Tools for representing the 3D space around a protein as a voxel grid of arbitrary resolution.

### Generating a csv file containing voxel centroids:

See below as executed in `if __name__ == __main__` block of collision_grid.py

```python
    import prody as pdy
    file_name = pdy.fetchPDB('1stp')
    mol = pdy.parsePDB(file_name)
    mol = mol.select('name CA')
    xs, ys, zs = [], [], []
    for atom in mol.getCoords():
        xs.append(atom[0])
        ys.append(atom[1])
        zs.append(atom[2])
    matrix = CollisionGrid(xs, ys, zs, 3, 10)
    matrix.write_voxel_coordinate_csv("test.csv", [(x,y,z) for x,y,z in zip(xs, ys, zs)])
```

Those centroids can be drawn in pymol by running cubes.py

The cubes.py script should be run from inside pymol. 
Just drag and drop it into the command line section on a mac or file->open the script on windows.

