'''
Square and Tetrahedra representations - Extended by Benjamin Fry
(bfry@g.harvard.edu). 

Should be run from inside pymol. Just file->open and select the script.

(c) 2013 Thomas Holder

License: BSD-2-Clause
'''
from pymol import cmd, cgo

def cgo_cube(x, y, z, r, color=None):
    # r *= 3 ** -.5
    color_line = [cgo.COLOR, *color] if color else None
    extend = [
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 0., 1.,
        cgo.VERTEX, x + r, y + r, z + r,
        cgo.VERTEX, x + r, y - r, z + r,
        cgo.VERTEX, x - r, y + r, z + r,
        cgo.VERTEX, x - r, y - r, z + r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 1., 0., 0.,
        cgo.VERTEX, x + r, y - r, z - r,
        cgo.VERTEX, x + r, y + r, z - r,
        cgo.VERTEX, x + r, y - r, z + r,
        cgo.VERTEX, x + r, y + r, z + r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 1., 0.,
        cgo.VERTEX, x + r, y + r, z - r,
        cgo.VERTEX, x - r, y + r, z - r,
        cgo.VERTEX, x + r, y + r, z + r,
        cgo.VERTEX, x - r, y + r, z + r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., 0., -1.,
        cgo.VERTEX, x - r, y - r, z - r,
        cgo.VERTEX, x - r, y + r, z - r,
        cgo.VERTEX, x + r, y - r, z - r,
        cgo.VERTEX, x + r, y + r, z - r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, -1., 0., 0.,
        cgo.VERTEX, x - r, y + r, z + r,
        cgo.VERTEX, x - r, y - r, z + r,
        cgo.VERTEX, x - r, y + r, z - r,
        cgo.VERTEX, x - r, y - r, z - r,
        cgo.END,
        cgo.BEGIN, cgo.TRIANGLE_STRIP,
        cgo.NORMAL, 0., -1., 0.,
        cgo.VERTEX, x - r, y - r, z + r,
        cgo.VERTEX, x + r, y - r, z + r,
        cgo.VERTEX, x - r, y - r, z - r,
        cgo.VERTEX, x + r, y - r, z - r,
        cgo.END,
    ]
    if color_line: 
        obj = color_line
        obj.extend(extend)
    else:
        obj = extend
    return obj

def draw_cubes_from_csv(file_name, color, name):
    with open(file_name, "r") as f:
        cubes = []
        for i,j in enumerate(f.readlines()):
            if i == 0:
                resolution = float(j.strip())
            else:
                x,y,z = j.strip().split(",")
                cube = cgo_cube(float(x),float(y),float(z),resolution/2, color)
                cubes.extend(cube)
        cmd.load_cgo(cubes, name, 1)

# Color options, you can extend by selecting normalized RGB values.
GREEN = (0.1, 1.0, 0.0)
RED = (1.0, 0.1, 0.0)
ORANGE = (1.0, 0.647, 0)
PURPLE = (0.8, 0.4, 0.9)
BLUE = (0.1, 0.0, 1.0)

# Load test.csv generated from running collision_grid.py
draw_cubes_from_csv("test.csv", BLUE, "1stp_voxelized")

# Set voxel transparency to 0.5
cmd.set("cgo_transparency", 0.5)
