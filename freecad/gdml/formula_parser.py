import re
from collections import defaultdict

class IllegalElement(Exception):
    def __init__(self, message):
        super().__init__(message)

def parse_chemical_formula(formula):
    recongized_elements = [
        'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
        'Cs', 'Ba', 'La',
        'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
        'Fr', 'Ra', 'Ac',
        'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf'
    ]

    def multiply_group(group_counts, multiplier):
        """Helper function to multiply the counts in a group by a multiplier."""
        for element in group_counts:
            group_counts[element] *= multiplier
        return group_counts

    def merge_counts(main_counts, group_counts):
        """Helper function to merge element counts from a group into the main dictionary."""
        for element, count in group_counts.items():
            main_counts[element] += count

    def parse_group(formula_part):
        """Recursively parse the formula part and return the counts."""
        element_counts = defaultdict(int)
        # Regular expression to find elements and their counts
        pattern = r'([A-Z][a-z]?)(\d*)|(\()|(\))(\d*)'
        stack = []
        group_counts = defaultdict(int)

        matches = re.finditer(pattern, formula_part)
        for match in matches:
            element, count, open_paren, close_paren, group_multiplier = match.groups()

            if element:
                if element not in recongized_elements:
                    raise IllegalElement(f"{element} not a valid element symbol")
                count = int(count) if count else 1
                element_counts[element] += count
            elif open_paren:
                # Push the current state onto the stack and start a new group
                stack.append((element_counts, group_counts))
                element_counts = defaultdict(int)
            elif close_paren:
                # Close the current group and multiply by the multiplier
                group_multiplier = int(group_multiplier) if group_multiplier else 1
                group_counts = multiply_group(element_counts, group_multiplier)
                element_counts, previous_group_counts = stack.pop()
                merge_counts(previous_group_counts, group_counts)
                element_counts.update(previous_group_counts)

        return element_counts

    # Parse the entire formula
    element_counts = parse_group(formula)

    return dict(element_counts)

# Example usage:
if __name__ == "__main__":
    import sys
    # formula = "H2(SO4)2"
    formula = sys.argv[1]
    print(formula)
    try:
        result = parse_chemical_formula(formula)
        print(result)  # Output will be {'H': 2, 'S': 2, 'O': 8}
    except Exception as e:
        print(e)

