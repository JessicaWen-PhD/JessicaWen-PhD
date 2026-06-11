"""
This script is used to create the SCF and environ input files for a series of charges using the results of a relax calculation,
as well as the bash files to run the calculations.

Usage:
    python charge-relax-structure.py --relax_name <relax_name> --output_name <output_name> --charge_start <charge_start> --charge_end <charge_end> --charge_increment <charge_increment>

Arguments:
    --relax_name: The name of the relax calculation.
    --output_name: The name of the output files.
    --charge_start: The start charge of the range.
    --charge_end: The end charge of the range (using Numpy's arange function, so exclusive of end).
    --charge_increment: The increment to go up between charge_start and charge_end

Example:
    python charge-relax-structure.py --relax_name Pt-Zundel --output_name Pt-Zundel --charge_start -0.5 --charge_end 1.1 --charge_increment 0.1
"""

# import relevant packages
import argparse
import re
from pathlib import Path
import numpy as np
from ase.io import read
from typing import List
from ase import Atoms

# define arguments
parser = argparse.ArgumentParser()
parser.add_argument("--relax_name", type=str, required=True)
parser.add_argument("--output_name", type=str, required=True)
parser.add_argument("--charge_start", type=float, required=True)
parser.add_argument("--charge_end", type=float, required=False)
parser.add_argument("--charge_increment", type=float, required=False)
args = parser.parse_args()
charge_range = np.arange(args.charge_start, args.charge_end, args.charge_increment)

POSITIONS_PLACEHOLDER = "__ATOMIC_POSITIONS_PLACEHOLDER__\n"


def _read_lines(path: Path) -> List[str]:
    with path.open("r") as f:
        return f.readlines()


def _write_lines(path: Path, lines: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.writelines(lines)


def extract_input_template(relax_in_path: Path) -> List[str]:
    """
    Returns QE relax input lines with the ATOMIC_POSITIONS block replaced by a single placeholder line.
    """
    raw = _read_lines(relax_in_path)

    # Replace the whole ATOMIC_POSITIONS block with a placeholder.
    try:
        pos_start = next(i for i, l in enumerate(raw) if "ATOMIC_POSITIONS" in l)
    except StopIteration:
        # Some inputs might not have ATOMIC_POSITIONS; just append placeholder.
        raw = raw + [POSITIONS_PLACEHOLDER]
    else:
        raw = raw[:pos_start] + [POSITIONS_PLACEHOLDER] + raw[pos_start + 1 :]

    return raw


def _upsert_namelist_var(lines: List[str], namelist: str, key: str, value: str) -> List[str]:
    """
    Insert or update `key = value` within a QE namelist section (&CONTROL, &SYSTEM, ...).
    `value` should include any needed quotes.
    """
    out = list(lines)
    nml_pat = re.compile(rf"^\s*&{re.escape(namelist)}\b", re.IGNORECASE)
    end_pat = re.compile(r"^\s*/\s*$")
    key_pat = re.compile(rf"^\s*{re.escape(key)}\s*=", re.IGNORECASE)

    in_nml = False
    nml_start_idx = None
    for i, line in enumerate(out):
        if not in_nml and nml_pat.search(line):
            in_nml = True
            nml_start_idx = i
            continue
        if in_nml and end_pat.search(line):
            # Not found; insert before "/" line.
            out.insert(i, f"    {key} = {value}\n")
            return out
        if in_nml and key_pat.search(line):
            out[i] = f"    {key} = {value}\n"
            return out

    if nml_start_idx is None:
        # Namelist not present; prepend a minimal one.
        return [f"&{namelist}\n", f"    {key} = {value}\n", "/\n"] + out

    # Namelist start found but no "/" end found; append.
    out.append(f"    {key} = {value}\n")
    return out


def _remove_namelist(lines: List[str], namelist: str) -> List[str]:
    """
    Remove an entire QE namelist block like:
      &IONS
        ...
      /
    If the block is not present, returns lines unchanged.
    """
    out: List[str] = []
    nml_pat = re.compile(rf"^\s*&{re.escape(namelist)}\b", re.IGNORECASE)
    end_pat = re.compile(r"^\s*/\s*$")

    skip = False
    for line in lines:
        if not skip and nml_pat.search(line):
            skip = True
            continue
        if skip:
            if end_pat.search(line):
                skip = False
            continue
        out.append(line)
    return out


def build_scf_template(input_path: Path, charge: float) -> List[str]:
    tmpl = extract_input_template(input_path)
    print(f"[debug] build_scf_template: loaded {len(tmpl)} lines from {input_path}")

    # Force relax settings in &CONTROL
    tmpl = _upsert_namelist_var(tmpl, "CONTROL", "calculation", "'relax'")
    tmpl = _upsert_namelist_var(tmpl, "CONTROL", "restart_mode", "'from_scratch'")
    tmpl = _upsert_namelist_var(tmpl, "CONTROL", "disk_io", "'nowf'")

    # Set total charge in &SYSTEM
    tmpl = _upsert_namelist_var(tmpl, "SYSTEM", "tot_charge", f"{charge:g}")
    print(f"[debug] build_scf_template: set tot_charge={charge:g}")

    return tmpl


def extract_relaxed_positions(relax_out_path: Path) -> Atoms:
    """
    Parses a QE relax output file and returns the lines between
    'Begin final coordinates' and 'End final coordinates'.
    """
    print(f"[debug] extract_relaxed_positions: reading {relax_out_path}")
    relaxed_structure = read(relax_out_path)
    print(f"[debug] extract_relaxed_positions: found {len(relaxed_structure)} atoms")
    return relaxed_structure


def create_relax_input_files(relax_name: str, output_name: str, charge: float) -> None:
    relax_in_path = Path(f"{relax_name}.relax.in")
    relax_out_path = Path(f"{relax_name}-success/{relax_name}.relax.out")
    print(f"[debug] create_relax_input_files: relax_in={relax_in_path}, charge={charge:.2f}")

    if not relax_in_path.exists():
        raise FileNotFoundError(f"Missing relax input file: {relax_in_path}")

    tmpl = build_scf_template(relax_in_path, charge)

    relaxed_atoms = extract_relaxed_positions(relax_out_path)
    pos_lines = ["ATOMIC_POSITIONS {angstrom}\n"]
    for symbol, pos in zip(relaxed_atoms.get_chemical_symbols(), relaxed_atoms.positions):
        pos_lines.append(f"{symbol}  {pos[0]:.10f}  {pos[1]:.10f}  {pos[2]:.10f}\n")

    # Drop original position lines from template; inject relaxed ones at placeholder
    in_old_positions = False
    new_tmpl = []
    for l in tmpl:
        if l == POSITIONS_PLACEHOLDER:
            new_tmpl.extend(pos_lines)
            in_old_positions = True
            continue
        if in_old_positions:
            # Skip lines that look like position data (element + 3 floats)
            if re.match(r"^\s*\w+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+", l):
                continue
            in_old_positions = False
        new_tmpl.append(l)
    tmpl = new_tmpl

    name = f"{output_name}{round(charge * 10):d}"
    out_path = Path(f"{name}/{output_name}{round(charge * 10):d}.relax.in")
    tmpl = _upsert_namelist_var(tmpl, "CONTROL", "outdir", f"\"./\"")
    tmpl = _upsert_namelist_var(tmpl, "CONTROL", "prefix", f"\"{name}\"")

    _write_lines(out_path, tmpl)
    print(f"[debug] create_relax_input_files: written to {out_path}")


_ENVIRON_TEMPLATE = """\
&ENVIRON
 verbose = 1
 environ_thr = 1.d-3
 environ_type = 'water'
 environ_restart = .FALSE.
 env_electrostatic = .TRUE.
/
&BOUNDARY
 solvent_mode = 'full'
 corespread(2) = 0.0
 solvent_radius = 2.6
 filling_threshold = 0.7
/
&ELECTROSTATIC
 pbc_correction = 'parabolic'
 pbc_dim = 2
 pbc_axis = 3
 tol = 1.d-13
/

EXTERNAL_CHARGES (angstrom)
{half_charge:.3f} 0. 0. {top_plane:.10f} 1.0 2 3
{half_charge:.3f} 0. 0. {bottom_plane:.10f} 1.0 2 3
"""


def create_environ_input_files(relax_name: str, output_name: str, charge: float) -> None:
    """
    Write an environ input file with external charge planes placed at max_z + 4
    and min_z - 4 of the relaxed atomic positions. External charge is half the
    total charge (split symmetrically above and below the slab).
    """

    relax_out_path = Path(f"{relax_name}-success/{relax_name}.relax.out")
    print(f"[debug] create_environ_input_files: relax_out={relax_out_path}, charge={charge:.2f}")

    relaxed_atoms = extract_relaxed_positions(relax_out_path)
    half_charge = -0.5 * charge
    print(f"[debug] create_environ_input_files: half_charge={half_charge:.4f}")

    z_coords = relaxed_atoms.positions[:, 2].tolist()

    if not z_coords:
        raise ValueError(f"Could not parse z-coordinates for {relax_name}")

    top_plane = max(z_coords) + 4.0
    bottom_plane = min(z_coords) - 4.0
    print(f"[debug] create_environ_input_files: z range [{min(z_coords):.4f}, {max(z_coords):.4f}]")
    print(f"[debug] create_environ_input_files: top_plane={top_plane:.4f}, bottom_plane={bottom_plane:.4f}")

    content = _ENVIRON_TEMPLATE.format(top_plane=top_plane, bottom_plane=bottom_plane, half_charge=half_charge)
    name = f"{output_name}{round(charge * 10):d}"
    out_path = Path(f"{name}/environ.in")
    _write_lines(out_path, [content])
    print(f"[debug] create_environ_input_files: written to {out_path}")


_BASH_TEMPLATE = """\
#!/bin/bash
#SBATCH -p cpuonly
#SBATCH -N 1
#SBATCH --ntasks-per-node=96
#SBATCH --job-name={name}
#SBATCH -t 47:30:00
#SBATCH -o {name}.out
#SBATCH -e {name}.err

ulimit -s unlimited
ulimit -a  # Print all limits for debugging
export OMP_NUM_THREADS=1

SECONDS=0
module purge
module load psc.allocations.user/1.0
module load intel-oneapi-compilers/2022.1.0 intel-oneapi-mkl/2022.1.0 intel-oneapi-mpi/2021.6.0

PW=/trace/group/dabo/shared/software/qe/qe-7.4.1/build/bin/pw.x
BASENAME="{name}"

mkdir -p $BASENAME
mpirun -np "${{SLURM_NTASKS}}" $PW --environ -in $BASENAME.relax.in > $BASENAME/$BASENAME.relax.out 2>&1

duration=$SECONDS
time=`date +%Y%m%d-%H%M%S`
echo "Relax is completed in $((duration / 60)) minutes and $((duration % 60)) seconds, date:${{time}}." >> $BASENAME/$BASENAME.relax.out
"""


def create_bash_files(output_name: str, charge: float) -> None:
    """
    Write a bash file that is able to run the input files with environ in QE.
    """
    name = f"{output_name}{round(charge * 10):d}"
    print(f"[debug] create_bash_files: generating script for job '{name}'")
    content = _BASH_TEMPLATE.format(name=name)
    out_path = Path(f"{name}/run_scf.sh")
    _write_lines(out_path, [content])
    print(f"[debug] create_bash_files: written to {out_path}")


if __name__ == "__main__":
    print(f"[debug] charge_range: {charge_range}")
    for charge in charge_range:
        print(f"\n[debug] ---- Processing charge={charge:.2f} ----")
        create_relax_input_files(args.relax_name, args.output_name, charge)
        create_environ_input_files(args.relax_name, args.output_name, charge)
        create_bash_files(args.output_name, charge)
