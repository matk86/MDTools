import numpy as np


from pymatgen import Molecule
from pymatgen.core.operations import SymmOp
from pymatgen.util.coord_utils import get_angle


__author__ = "Kiran Mathew, Brandon Wood"


class Polymer(object):
    def __init__(self, start_monomer, s_head, s_tail,
                 monomer, head, tail,
                 end_monomer, e_head, e_tail,
                 n_units, link_distance=1.0, linear_chain=False):
        """
        Args:
            start_monomer (Molecule): Starting molecule
            s_head (int): starting atom index of the start_monomer molecule
            s_tail (int): tail atom index of the start_monomer
            monomer (Molecule): The monomer
            head (int): index of the atom in the monomer that forms the head
            tail (int): tail atom index. monomers will be connected from
                tail to head
            end_monomer (Molecule): Terminal molecule
            e_head (int): starting atom index of the end_monomer molecule
            e_tail (int): tail atom index of the end_monomer
            n_units (int): number of monomer units excluding the start and
                terminal molecules
            link_distance (float): distance between consecutive monomers
            linear_chain (bool): linear or random walk polymer chain
        """
        self.start = s_head
        self.end = s_tail
        self.monomer = monomer
        self.n_units = n_units
        self.link_distance = link_distance
        self.linear_chain = linear_chain
        # translate monomers so that head atom is at the origin
        start_monomer.translate_sites(range(len(start_monomer)),
                                      - monomer.cart_coords[s_head])
        monomer.translate_sites(range(len(monomer)),
                                - monomer.cart_coords[head])
        end_monomer.translate_sites(range(len(end_monomer)),
                                    - monomer.cart_coords[e_head])
        self.mon_vector = monomer.cart_coords[tail] - monomer.cart_coords[head]
        self.moves = {1: [1, 0, 0], 2: [0, 1, 0], 3: [0, 0, 1], 4: [-1, 0, 0],
                      5: [0, -1, 0], 6: [0, 0, -1]}
        self.prev_move = 1
        # places the start monomer at the beginning of the chain
        self.molecule = start_monomer.copy()
        self.length = 0
        # create the chain
        self._create(self.monomer, self.mon_vector)
        # terminate the chain with the end_monomer
        self.n_units += 1
        end_mon_vector = end_monomer.cart_coords[e_tail] - \
                         end_monomer.cart_coords[e_head]
        self._create(end_monomer, end_mon_vector)

    def _create(self, monomer, mon_vector):
        """
        create the polymer from the monomer

        Args:
            monomer (Molecule)
            mon_vector (numpy.array): molecule vector that starts from the
                start atom index to the end atom index
        """
        while self.length != self.n_units:
            if self.linear_chain:
                move_direction = np.array(mon_vector) / np.linalg.norm(
                    mon_vector)
            else:
                move_direction = self._next_move_direction()
            self._add_monomer(monomer.copy(), mon_vector, move_direction)

    def _next_move_direction(self):
        """
        pick a move at random from the list of moves
        """
        move = np.random.randint(1, 7)
        while self.prev_move == (move + 3) % 6:
            move = np.random.randint(1, 7)
        self.prev_move = move
        #print "move", move
        return np.array(self.moves[move])

    def _align_monomer(self, monomer, mon_vector, move_direction):
        """
        rotate the monomer so that it is aligned along the move direction

        Args:
            monomer (Molecule)
            mon_vector (numpy.array): molecule vector that starts from the
                start atom index to the end atom index
            move_direction (numpy.array): the direction of the polymer chain
                extension
        """
        axis = np.cross(mon_vector, move_direction)
        origin = monomer[self.start].coords
        angle = get_angle(mon_vector, move_direction)
        #print "axis, origin, angle", axis, origin, angle
        op = SymmOp.from_origin_axis_angle(origin, axis, angle)
        monomer.apply_operation(op)

    def _add_monomer(self, monomer, mon_vector, move_direction):
        """
        extend the polymer molecule by adding a monomer along mon_vector direction

        Args:
            monomer (Molecule): monomer molecule
            mon_vector (numpy.array): monomer vector that points from head to tail.
            move_direction (numpy.array): direction along which the monomer
                will be positioned
        """
        translate_by = self.molecule.cart_coords[
                           self.end] + self.link_distance * move_direction
        monomer.translate_sites(range(len(monomer)), translate_by)
        if not self.linear_chain:
            self._align_monomer(monomer, mon_vector, move_direction)
        # add monomer if there are no crossings
        does_cross = False
        for i, site in enumerate(monomer):
            try:
                self.molecule.append(site.specie, site.coords,
                                     properties=site.properties)
            except:
                does_cross = True
                polymer_length = len(self.molecule)
                self.molecule.remove_sites(
                    range(polymer_length - i, polymer_length))
                break
        if not does_cross:
            self.length += 1
            self.end += len(self.monomer)
            #print "length: ", self.length
