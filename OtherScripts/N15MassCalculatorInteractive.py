"""
Description: Script for calculating the N14 and N15 mass from a sequence
"""
import pyteomics
from pyteomics import mass


def calculate_mass(seq: str, ch: int) -> tuple:
    # Get the composition
    compositions = list(pyteomics.mass.isotopologues(seq, elements_with_isotopes='N'))
    n14_comp = compositions[0]
    n15_comp = compositions[-1]
    # Get the m/z
    n14_mass = round(pyteomics.mass.calculate_mass(n14_comp, ch), 3)
    n15_mass = round(pyteomics.mass.calculate_mass(n15_comp, ch), 3)
    return n14_mass, n15_mass


if __name__ == '__main__':
    print("Welcome to the N14 and N15 mass calculator")
    new_sequence: bool = True
    while new_sequence:
        sequence = input("Enter sequence: ").upper()
        charge = int(input("Enter charge: "))
        (n14, n15) = calculate_mass(sequence, charge)

        print(f"For sequence {sequence} ({charge}+) the N14 m/z is {n14} and N15 m/z is {n15}")
        new_sequence = input("Enter a new sequence (y/n): ").lower() == 'y'
    print('Bye bye')
