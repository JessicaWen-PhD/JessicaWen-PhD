"""
This script is used to create the SCF and environ input files for a NEB calculation,
as well as the bash files to run the calculations.

Usage:
    python charge-neb-images.py --neb_name <neb_name> --output_name <output_name> --charge <charge>

Arguments:
    --neb_name: The name of the NEB calculation.
    --output_name: The name of the output files.
    --charge: The charge of the system.

Example:
    python charge-neb-images.py --neb_name Tafel-one-proton-charged --output_name Tafel-1proton-1 --charge 0.1
"""

# import relevant packages
import argparse
import re
import shutil
from pathlib import Path
from typing import List

# define arguments
parser = argparse.ArgumentParser()
parser.add_argument("--neb_name", type=str, required=True)
parser.add_argument("--output_name", type=str, required=True)
parser.add_argument("--charge", type=float, required=True)
args = parser.parse_args()

POSITIONS_PLACEHOLDER = "__ATOMIC_POSITIONS_PLACEHOLDER__\n"


def _read_lines(path: Path) -> List[str]:
    with path.open("r") as f:
        return f.readlines()


def _write_lines(path: Path, lines: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.writelines(lines)


def extract_number_of_images(neb_in_path: Path) -> int:
    for line in _read_lines(neb_in_path):
        if "num_of_images" in line:
            m = re.search(r"num_of_images\s*=\s*([0-9]+)", line)
            if m:
                return int(m.group(1))
            # fallback: "num_of_images     = 15,"
            return int(re.sub(r"[^0-9]", "", line.split("=")[-1]))
    raise ValueError(f"Could not find num_of_images in {neb_in_path}")


def extract_nat(neb_in_path: Path) -> int:
    # nat is in &SYSTEM namelist: e.g. "nat              = 16"
    for line in _read_lines(neb_in_path):
        if re.search(r"\bnat\b", line) and "=" in line:
            m = re.search(r"\bnat\s*=\s*([0-9]+)", line)
            if m:
                return int(m.group(1))
    raise ValueError(f"Could not find nat in {neb_in_path}")


def extract_engine_input_template(neb_in_path: Path) -> List[str]:
    """
    Returns QE engine input lines (between BEGIN_ENGINE_INPUT and END_ENGINE_INPUT),
    with the NEB BEGIN_POSITIONS block replaced by a single placeholder line.
    """
    raw = _read_lines(neb_in_path)

    try:
        start = next(i for i, l in enumerate(raw) if "BEGIN_ENGINE_INPUT" in l) + 1
        end = next(i for i, l in enumerate(raw) if "END_ENGINE_INPUT" in l)
    except StopIteration as e:
        raise ValueError(f"Could not locate engine input block in {neb_in_path}") from e

    engine = raw[start:end]

    # Replace the whole BEGIN_POSITIONS ... END_POSITIONS block with a placeholder.
    try:
        pos_start = next(i for i, l in enumerate(engine) if "BEGIN_POSITIONS" in l)
        pos_end = next(i for i, l in enumerate(engine) if "END_POSITIONS" in l)
    except StopIteration:
        # Some inputs might not have BEGIN_POSITIONS; just append placeholder.
        engine = engine + [POSITIONS_PLACEHOLDER]
    else:
        engine = engine[:pos_start] + [POSITIONS_PLACEHOLDER] + engine[pos_end + 1 :]

    return engine


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


def build_scf_template(neb_in_path: Path, charge: float) -> List[str]:
    tmpl = extract_engine_input_template(neb_in_path)

    # Force SCF settings in &CONTROL
    tmpl = _upsert_namelist_var(tmpl, "CONTROL", "calculation", "'scf'")
    tmpl = _upsert_namelist_var(tmpl, "CONTROL", "restart_mode", "'from_scratch'")
    tmpl = _upsert_namelist_var(tmpl, "CONTROL", "disk_io", "'nowf'")

    # Set total charge in &SYSTEM
    tmpl = _upsert_namelist_var(tmpl, "SYSTEM", "tot_charge", f"{charge:g}")

    # SCF: remove ionic dynamics namelist entirely
    tmpl = _remove_namelist(tmpl, "IONS")
    return tmpl


def extract_image_positions_from_crd(crd_path: Path, nat: int) -> List[List[str]]:
    """
    Parses a QE NEB .crd file formatted like:
      FIRST_IMAGE / INTERMEDIATE_IMAGE / LAST_IMAGE
      ATOMIC_POSITIONS (angstrom)
        <nat lines>
    Returns a list of per-image blocks, each a list[str] including the ATOMIC_POSITIONS line.
    """
    lines = _read_lines(crd_path)
    image_header_pat = re.compile(r"^\s*(FIRST_IMAGE|INTERMEDIATE_IMAGE|LAST_IMAGE)\s*$")
    atomic_pos_pat = re.compile(r"^\s*ATOMIC_POSITIONS\b", re.IGNORECASE)

    blocks: List[List[str]] = []
    i = 0
    while i < len(lines):
        if not image_header_pat.match(lines[i]):
            i += 1
            continue

        # Expect ATOMIC_POSITIONS next (allow blank lines in between)
        j = i + 1
        while j < len(lines) and lines[j].strip() == "":
            j += 1
        if j >= len(lines) or not atomic_pos_pat.match(lines[j]):
            raise ValueError(f"Unexpected .crd format near line {i+1} in {crd_path}")

        block = [lines[j]]
        start_coords = j + 1
        coords = lines[start_coords : start_coords + nat]
        if len(coords) != nat:
            raise ValueError(f"Not enough coordinate lines for nat={nat} near line {j+1} in {crd_path}")
        block.extend(coords)
        blocks.append(block)

        i = start_coords + nat

    if not blocks:
        raise ValueError(f"Could not find any image blocks in {crd_path}")
    return blocks

def create_scf_input_files(neb_name: str, output_name: str, charge: float) -> None:
    neb_in_path = Path(f"{neb_name}.neb.in")
    crd_path = Path(f"{neb_name}.crd")

    if not neb_in_path.exists():
        raise FileNotFoundError(f"Missing NEB input file: {neb_in_path}")
    if not crd_path.exists():
        raise FileNotFoundError(f"Missing coordinates file: {crd_path}")

    nat = extract_nat(neb_in_path)
    declared_images = extract_number_of_images(neb_in_path)
    image_blocks = extract_image_positions_from_crd(crd_path, nat)

    if declared_images != len(image_blocks):
        raise ValueError(
            f"num_of_images in {neb_in_path} is {declared_images}, but {crd_path} contains {len(image_blocks)} image blocks"
        )

    tmpl = build_scf_template(neb_in_path, charge)

    for idx, positions_block in enumerate(image_blocks, start=1):
        out_path = Path(f"{output_name}-{idx}/{output_name}-{idx}.scf.in")
        rendered: List[str] = []
        # Change outdir and prefix to the output name
        image_tmpl = _upsert_namelist_var(tmpl, "CONTROL", "outdir", f"\"./\"")
        image_tmpl = _upsert_namelist_var(image_tmpl, "CONTROL", "prefix", f"\"{output_name}-{idx}\"")

        for line in image_tmpl:
            if line == POSITIONS_PLACEHOLDER:
                rendered.extend(positions_block)
            else:
                rendered.append(line)

        _write_lines(out_path, rendered)

_ENVIRON_TEMPLATE = """\
&ENVIRON
 verbose = 2
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
 pbc_correction = 'none'
 pbc_dim = 2
 pbc_axis = 3
 tol = 1.d-13
/

EXTERNAL_CHARGES (angstrom)
{half_charge:.3f} 0. 0. {top_plane:.10f} 1.0 2 3
{half_charge:.3f} 0. 0. {bottom_plane:.10f} 1.0 2 3
"""


def create_environ_input_files(neb_name: str, output_name: str, charge: float) -> None:
    """
    For each NEB image, write an environ input file with external charge planes
    placed at max_z + 4 and min_z - 4 of the image's atomic positions. External
    charge is given by half of the charge in the QE input file.
    """
    neb_in_path = Path(f"{neb_name}.neb.in")
    crd_path = Path(f"{neb_name}.crd")
    half_charge = -0.5 * charge

    if not neb_in_path.exists():
        raise FileNotFoundError(f"Missing NEB input file: {neb_in_path}")
    if not crd_path.exists():
        raise FileNotFoundError(f"Missing coordinates file: {crd_path}")

    nat = extract_nat(neb_in_path)
    image_blocks = extract_image_positions_from_crd(crd_path, nat)

    for idx, positions_block in enumerate(image_blocks, start=1):
        # positions_block[0] is the ATOMIC_POSITIONS header; remaining lines are coords
        z_coords = []
        for line in positions_block[1:]:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    z_coords.append(float(parts[3]))
                except ValueError:
                    pass

        if not z_coords:
            raise ValueError(f"Could not parse z-coordinates for image {idx}")

        top_plane = max(z_coords) + 4.0
        bottom_plane = min(z_coords) - 4.0

        content = _ENVIRON_TEMPLATE.format(top_plane=top_plane, bottom_plane=bottom_plane, half_charge=half_charge)
        out_path = Path(f"{output_name}-{idx}/environ.in")
        _write_lines(out_path, [content])

_BASH_TEMPLATE = """\
#!/bin/bash
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --job-name={name}
#SBATCH --account=chm230047
#SBATCH --time=2:00:00

cd $SLURM_SUBMIT_DIR

source /anvil/projects/x-chm230047/Shared/software/setup_intel_impi_env.sh

QE_DIR=/anvil/projects/x-chm230047/Shared/software/qe/qe-7.5-intel_impi/build/bin
BASENAME="{name}"

echo "Job started on $(hostname) at $(date)"
echo "Running with $SLURM_NTASKS MPI tasks"
echo

# === QE calculations ===

echo "Running QE: scf + sccs on $BASENAME"
mkdir -p $BASENAME
srun --mpi=pmi2 --kill-on-bad-exit=1 -n $SLURM_NTASKS $QE_DIR/pw.x --environ -in $BASENAME.scf.in > $BASENAME/$BASENAME.scf.out

duration=$SECONDS
time=`date +%Y%m%d-%H%M%S`
echo "SCF is completed in $((duration / 60)) minutes and $((duration % 60)) seconds, date:${{time}}." >> $BASENAME/$BASENAME.scf.out
"""

def create_bash_files(neb_name: str, output_name: str):
    """
    For each NEB image, write a bash file that is able to run the input files with environ in QE.
    """
    neb_in_path = Path(f"{neb_name}.neb.in")
    crd_path = Path(f"{neb_name}.crd")

    if not neb_in_path.exists():
        raise FileNotFoundError(f"Missing NEB input file: {neb_in_path}")
    if not crd_path.exists():
        raise FileNotFoundError(f"Missing coordinates file: {crd_path}")

    nat = extract_nat(neb_in_path)
    image_blocks = extract_image_positions_from_crd(crd_path, nat)

    for idx, positions_block in enumerate(image_blocks, start=1):
        content = _BASH_TEMPLATE.format(name=f"{output_name}-{idx}")
        out_path = Path(f"{output_name}-{idx}/run_scf.sh")
        _write_lines(out_path, [content])


if __name__ == "__main__":
    create_scf_input_files(args.neb_name, args.output_name, args.charge)
    create_environ_input_files(args.neb_name, args.output_name, args.charge)
    create_bash_files(args.neb_name, args.output_name)