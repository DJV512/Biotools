from dataclasses import dataclass
from . import constants as btc


class Protein:

    def __init__(self, name:str, sequence:str, allow_x=False):
        self.name = name
        self.sequence = sequence.upper()
        self._validate_protein(allow_x = allow_x)


    def __str__(self):
        return self.name


    def __eq__(self, other):
        if isinstance(other, Protein):
            return self.name == other.name and self.sequence == other.sequence 


    def _validate_protein(self, allow_x):
        '''
        Only the normal 20 amino acids are allowed
        
        If allow_x is True, 'X' is also allowed.
        
        'X' is used as an alignment anchor in Multiple Sequence Alignment.
        '''

        allowed_aas = set(btc.AMINO_ACIDS)
        if allow_x:
            allowed_aas.add("X")
        residues = set(self.sequence)
        if not residues.issubset(allowed_aas):
            raise ValueError(f"{self.sequence} is not valid a Protein sequence.")


    def charge(self) -> int:
        '''
        Calculates and returns its overall charge (at pH 7).
        '''

        charge = 0
        for aa in self.sequence:
            if aa in ["K", "R"]:
                charge += 1
            elif aa in ["D", "E"]:
                charge -= 1
        return charge


    def mw(self) -> float:
        '''
        Given a protein, calculates and returns its molecular weight.
        '''
        self._validate_protein(allow_x=False)
        weight = 18
        for aa in self.sequence:
            weight += btc.MW[aa]
        return weight


    def ec(self) -> int:
        '''
        Given a protein sequence, calculates and returns the extinction coefficient, in cm-1/M.
        '''

        ec = 0
        for aa in self.sequence:
            if aa == "W":
                ec += 5690
            elif aa == "Y":
                ec += 1280
            elif aa == "C":
                ec += 120
        return ec


    def glycosylation(self) -> list:
        '''
        Given a protein sequence, finds all potential N-linked glycosylations sites, denoted by the motif N-X-S/T, where X is any amino acid except proline.
        '''

        sites = []
        length = len(self.sequence)
        for i, aa in enumerate(self.sequence):
            if aa == "N" and i < length - 2:
                aa2 = self.sequence[i+1]
                if aa2 != "P":
                    aa3 = self.sequence[i+2]
                    if aa3 in ["T", "S"]:
                        sites.append(i)
        return sites


@dataclass
class AlignedProteinPair:
    protein1: Protein
    protein2: Protein
    align1: list[str]
    align2: list[str]
    distance: float


@dataclass
class Cluster:
    id: str
    name: str
    left: Protein | Cluster
    right: Protein | Cluster
    members: list[Protein]

    OFFSET = 8

    def __str__(self):
        members = ", ".join(member.name for member in self.members)
        return f"Cluster(name={self.name}, left={self.left}, right={self.right}, members={members})"


    def swap_nodes(self):
        left_node = self.left
        right_node = self.right
        self.left = right_node
        self.right = left_node


    def _max_height(self):
        '''
        Determines the number of levels on a guide tree.
        '''
        if isinstance(self.left, Cluster):
            left_height =  1 + self.left._max_height()
        else:
            left_height = 1

        if isinstance(self.right, Cluster):
            right_height = 1 + self.right._max_height()
        else:
            right_height = 1

        return max(left_height, right_height)


    def _find_coords(self):
        '''
        Finds the coordinates of Clusters for the purposes of printing the guide tree.
        '''
        coords = {}
        coords = self._find_x_coords(coords)
        coords = self._find_y_coords(coords)
        coords = self._find_parent_coords(coords)
        return coords
        

    def _find_x_coords(self, coords):
        '''
        Finds Cluster x coordinates.
        '''
        x = 0
        for member in self.members:
            coords[member.name] = (x, 0)
            x += Cluster.OFFSET
        return coords


    def _find_y_coords(self, coords, depth = OFFSET):
        '''
        Finds Cluster y coordinates.
        '''
        if isinstance(self.left, Protein):
            x, _ = coords[self.left.name]
            coords[self.left.name] = (x, depth)
        else:
            self.left._find_y_coords(coords, depth + Cluster.OFFSET)

        if isinstance(self.right, Protein):
            x, _ = coords[self.right.name]
            coords[self.right.name] = (x, depth)
        else:
            self.right._find_y_coords(coords, depth + Cluster.OFFSET)
        
        return coords


    def _find_parent_coords(self, coords):
        '''
        Find coordinates of Cluster parents.
        '''
        if isinstance(self.left, Cluster):
            coords = self.left._find_parent_coords(coords)

        if isinstance(self.right, Cluster):
            coords = self.right._find_parent_coords(coords)

        x1, y = coords[self.left.name]
        x2, _ = coords[self.right.name]
        new_x = (x1 + x2) / 2
        new_y = y - Cluster.OFFSET
        coords[f"{self.left.name}-{self.right.name}"] = (new_x, new_y)

        return coords
        
    
    def print_tree(self):
        '''
        Method to pretty print a Cluster Tree for easy viewing.
        '''
        left_x_offset = 10
        right_x_offset = 20
        coords = self._find_coords()
        max_x = 0
        max_y = 0
        new_coords = {}
        for name, (x, y) in coords.items():
            if "-" in name and y == 0:
                new_name = "root"
            elif "-" in name:
                new_name = "*"
            else:
                new_name = name

            new_x = round(x)
            new_y = round(y)
            max_x = max(round(x), max_x)
            max_y = max(round(y), max_y)

            new_coords[(new_y, new_x + left_x_offset)] = new_name
        
        grid = [[" " for _ in range(max_x + left_x_offset + right_x_offset)] for _ in range(max_y+1)]

        for y in range(max_y + 1):
            for x in range(max_x + left_x_offset + 20):
                if (y, x) in new_coords:
                    length = len(new_coords[(y,x)])
                    start = x - length // 2
                    for i, char in enumerate(new_coords[(y, x)]):
                        grid[y][start + i + left_x_offset] = char
                elif (y - 1, x + 1) in new_coords and new_coords[(y-1, x + 1)] in ["*", "/", "root"] and y % Cluster.OFFSET != 0:
                    grid[y][x + left_x_offset] = "/"
                    new_coords[(y, x)] = "/"
                elif (y-1, x-1) in new_coords and new_coords[(y-1, x-1)] in ["*", "\\", "root"] and y % Cluster.OFFSET != 0:
                    grid[y][x + left_x_offset] = "\\"
                    new_coords[(y, x)] = "\\"

        for y in grid:
            print("".join(y))


@dataclass
class Profile:
    name: list[str]
    alignment: list[list[str]]

    def __str__(self):
        text = ""
        for name, align in zip(self.name, self.alignment):
            text += f"{name}, {align}\n"

        return text


@dataclass
class ParsedProfile:
    name: list[str]
    alignment: list[list[tuple[str, str]]]