#!/usr/bin/env python3
import numpy as np
#  from pyrosetta import *


class CollisionGrid:
    """Collision Grid data structure. Represents the volume around a pose as collection of voxels. Can be used to
    check if structures are colliding in linear time after initialization. """

    def __init__(self, xs=None, ys=None, zs=None, resolution=None, buffer_distance=None, class_instance=None):
        if class_instance == None:
            self.default_constructor(xs, ys, zs, resolution, buffer_distance)
        else:
            self.copy_constructor(class_instance)

    def default_constructor(self, xs, ys, zs, resolution, buffer_distance):
        """
        Default constructor to be used to create CollisionGrid objects.
        Pass a list of x,y,z coordinates corresponding to the largest object you want to voxelize.
        For example the xyz coordinates of all the atoms in a protein.
        """
        # Throw an error if any of the parameters are not initialized.
        if resolution == None or buffer_distance == None:
            raise ValueError("Initialize CollisionGrid xs, ys, zs, resolution, and buffer_distance")

        self.resolution = resolution
        self.buffer_distance = buffer_distance
        self.num_cells_populated = 0

        # Generate the (discretized) ranges of coordinates that the pose can occupy
        self.x_range = np.arange(int(np.floor(min(xs) - buffer_distance)),
                                 int(np.ceil(max(xs) + buffer_distance + resolution)), resolution)
        self.y_range = np.arange(int(np.floor(min(ys) - buffer_distance)),
                                 int(np.ceil(max(ys) + buffer_distance + resolution)), resolution)
        self.z_range = np.arange(int(np.floor(min(zs) - buffer_distance)),
                                 int(np.ceil(max(zs) + buffer_distance + resolution)), resolution)

        # Calculate the difference between largest and smallest coordinates
        # Add a buffer to allow space for glycans.
        deltx = (self.x_range[-1] - self.x_range[0])
        delty = (self.y_range[-1] - self.y_range[0])
        deltz = (self.z_range[-1] - self.z_range[0])

        # Initialize the matrix with all zeros.
        self.space_matrix = np.zeros((int(deltx / resolution), int(delty / resolution), int(deltz / resolution)))
        print(f"Built Matrix of Shape {self.space_matrix.shape}, at resolution {self.resolution} Angstrom(s),"
              f" with buffer of {self.buffer_distance} Angstrom(s) in all directions")

    def copy_constructor(self, class_instance):
        """Copy constructor used to create deep-copied CollisionGrid objects
        from existing ones."""
        self.resolution = class_instance.resolution
        self.buffer_distance = class_instance.buffer_distance
        self.num_cells_populated = class_instance.num_cells_populated
        self.x_range = class_instance.x_range.copy()
        self.y_range = class_instance.y_range.copy()
        self.z_range = class_instance.z_range.copy()
        self.space_matrix = class_instance.space_matrix.copy()

    @classmethod
    def copy_from_existing(cls, class_to_copy):
        """Class method that calls the copy constructor on the class passed in.
        Run with CollisionGrid.copy_from_existing(class_to_copy). Returns new
        copied class."""
        copy = cls(class_instance=class_to_copy)
        return copy

    def calc_matrix_idcs_from_coord(self, coordinate):
        """Given a list of 3D coordinates, returns the x, y, and z index of the
        space matrix that contains that point in discretized space"""
        x, y, z = coordinate
        xidx = int(np.floor((x - self.x_range[0]) / self.resolution))
        yidx = int(np.floor((y - self.y_range[0]) / self.resolution))
        zidx = int(np.floor((z - self.z_range[0]) / self.resolution))
        return xidx, yidx, zidx

    def query_matrix(self, list_of_coordinates):
        """Counts number of unique, occupied boxes that each coordinate falls in.
        Returns count."""
        atom_query_result = 0
        all_idcs = set()  # Create a set to hold unique indices
        for coordinate in list_of_coordinates:
            xidx, yidx, zidx = self.calc_matrix_idcs_from_coord(coordinate)
            all_idcs.add((xidx, yidx, zidx))  # Add each index, set ignores duplicates
        for u_xidx, u_yidx, u_zidx in list(all_idcs):  # Loop through unique indices and add result
            atom_query_result += self.space_matrix[u_xidx, u_yidx, u_zidx]
        return atom_query_result

    def populate_matrix(self, list_of_coordinates):
        """Populates matrix by marking every voxel a member of list of
        coordinates falls within as occupied (1). Returns number of voxels
        populated and updates num_cells_populated."""
        num_populated = 0  # Keep track of number of unique coordinates populated each call
        for coordinate in list_of_coordinates:
            xidx, yidx, zidx = self.calc_matrix_idcs_from_coord(coordinate)
            try:
                if self.space_matrix[xidx, yidx, zidx] != 1:
                    self.space_matrix[xidx, yidx, zidx] = 1
                    num_populated += 1  # If mtx isn't filled at position, fill and increment num_populated
            except IndexError:
                print("Warning! Tried to populate a cell out of bounds. This is probably okay but you may need to "
                      "increase the buffer size if occurring during glycan sampling step.")

        self.num_cells_populated += num_populated  # Update class variable.
        return num_populated

    def query_distance_around_coordinate(self, coord, distance):
        """Given a distance to check around a coordinate, returns number of
        filled cells within distance of the coordinate"""
        num_collisions = 0
        num_cells = int(np.ceil(distance / self.resolution))
        xidx, yidx, zidx = self.calc_matrix_idcs_from_coord(coord)
        for i in range(-num_cells, num_cells + 1):
            for j in range(-num_cells, num_cells + 1):
                for k in range(-num_cells, num_cells + 1):
                    try:
                        x, y, z = self.get_voxel_centers_from_list([[coord[0] + i, coord[1] + j, coord[2] + k]])[0]
                        # Only check if center of voxel is within distance from point of interest
                        sq_distance_btwn = (coord[0] - x) ** 2 + (coord[1] - y) ** 2 + (coord[2] - z) ** 2
                        if sq_distance_btwn <= distance ** 2:
                            num_collisions += self.space_matrix[xidx + i, yidx + j, zidx + k]
                    except IndexError:
                        pass
        return num_collisions

    def get_voxel_centers_from_list(self, list_of_coordinates):
        """Calculates the center of the cube representing the range of
         each index in the space matrix. Returns a list of unique tuples."""

        output_set = set()  # Use a set to only get unique centers.
        for coordinate in list_of_coordinates:
            # Calculate the index in range array for a given point.
            xidx, yidx, zidx = self.calc_matrix_idcs_from_coord(coordinate)

            # Average the position of index in the range array with that
            # index + resolution to calculate center of each range
            xctr = (self.x_range[xidx] + (self.x_range[xidx] + self.resolution)) / 2
            yctr = (self.y_range[yidx] + (self.y_range[yidx] + self.resolution)) / 2
            zctr = (self.z_range[zidx] + (self.z_range[zidx] + self.resolution)) / 2
            posn = (xctr, yctr, zctr)
            output_set.add(posn)  # insert center to output_set if unique
        return list(output_set)

    def write_voxel_coordinate_csv(self, file_name, list_of_coordinates):
        """Writes the output of get_voxel_centers_from_list method to a csv
        file (file_name). Prefixes csv with the resolution for plotting with
        cubes.py helper script."""
        with open(file_name, "w") as f:
            f.write(f"{self.resolution}\n")
            for x, y, z in self.get_voxel_centers_from_list(list_of_coordinates):
                f.write(f"{x},{y},{z}\n")

if __name__ == "__main__":

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


