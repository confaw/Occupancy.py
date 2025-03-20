import numpy as np
import re
from pymatgen.io.vasp import Vasprun, Outcar, Potcar, Poscar, Chgcar
from pymatgen.electronic_structure.core import Spin

outcar = Outcar("OUTCAR")
fermi = outcar.efermi
print("Fermi energy from OUTCAR:", fermi)
vasprun = Vasprun("vasprun.xml", parse_eigen=True)
bs = vasprun.get_band_structure(line_mode=False)
bs.efermi = fermi

first_kpt_idx = 0
occupied_first = np.sum(bs.bands[Spin.up][:, first_kpt_idx] < fermi)
print("Number of occupied bands at the first k-point (Spin.up):", occupied_first)

last_kpt_idx = -1
occupied_last = np.sum(bs.bands[Spin.up][:, last_kpt_idx] < fermi)
print("Number of occupied bands at the last k-point (Spin.up):", occupied_last)

potcar = Potcar.from_file("POTCAR")
poscar = Poscar.from_file("POSCAR")
structure = poscar.structure
composition = structure.composition

elements_in_poscar_order = [str(sp) for sp in structure.composition.elements]
total_valence_electrons = 0

for element, potcar_entry in zip(elements_in_poscar_order, potcar):
    zval = None
    lines = potcar_entry.data.split("\n")
    for line in lines:
        match = re.search(r"ZVAL\s*=\s*([\d\.]+)", line)
        if match:
            zval = float(match.group(1))
            break

    if zval is None:
        raise ValueError(f"Could not parse ZVAL from the POTCAR entry for element {element}.")
    count = composition[element]
    valence_contribution = zval * count
    total_valence_electrons += valence_contribution

    print(f"Element: {element}, ZVAL: {zval}, Count: {count}, "
          f"Valence Contribution: {valence_contribution}")

half_valence_electrons = total_valence_electrons/2

print(f"Total valence electrons in the structure: {total_valence_electrons}")
print(f"1/2 of the total valence electrons in the structure: {half_valence_electrons}")

def parse_outcar_magnetization(outcar_file="OUTCAR"):
    """
    Manually parse the final magnetization from the OUTCAR by searching 
    for the line containing 'number of electron' and 'magnetization'.
    Returns (number_of_electrons, magnetization) if found, or (None, None) otherwise.
    """
    number_of_electrons = None
    magnetization = None

    # This regex will look for something like:
    # number of electron\s+(\S+)\s+magnetization\s+(\S+)
    # capturing the numeric values.
    pattern = re.compile(
        r"number\s+of\s+electron\s+([\d\.\-Ee]+)\s+magnetization\s+([\d\.\-Ee]+)",
        re.IGNORECASE
    )

    with open(outcar_file, "r") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                number_of_electrons = float(match.group(1))
                magnetization = float(match.group(2))
                # Because VASP may print multiple iterations, 
                # we keep reading until the end to capture the *last* occurrence.

    return number_of_electrons, magnetization

if __name__ == "__main__":
    ne, mag = parse_outcar_magnetization("OUTCAR")
    if mag is not None:
        print(f"Found final number of electrons = {ne:.6f}, magnetization = {mag:.6f} µB")
    else:
        print("Could not find magnetization line in OUTCAR.")
