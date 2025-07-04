# RepairPBC - GROMACS Trajectory PBC Repair Tool

A command-line tool to repair periodic boundary conditions (PBC) in GROMACS trajectory files (.xtc) using their associated topology files (.gro).

## What it does
GROMACS simulations typically employ periodic boundary conditions (PBC) to constrain simulation box size. This means that if atom positions are extracted from the trajectories directly, proteins may appear "broken" if they cross outside the box boundary and "jump" across to the other side of the simulation box. 
This tool analyzes GROMACS trajectory files and corrects atoms that have "jumped" across periodic boundaries during molecular dynamics simulations. It uses an efficient array-based algorithm to identify and correct these jumps, ensuring your protein structure remains continuous throughout the trajectory.

## Important Limitation

This tool only works with single, contiguous protein chains. It will automatically validate your input files and refuse to process:

- Multiple protein chains (different molecules)
- Non-contiguous protein segments (gaps in residue numbering)
- Complex protein systems with multiple chains or domains

The algorithm assumes atoms are sequentially ordered in a single polypeptide chain. Using it on multi-chain systems will either throw an error or produce incorrect results.

## Installation

### Option 1: Pre-built Executable (Recommended for End Users)
For users who want to run the tool without installing Python or dependencies:

1. Go to the [Releases page](https://github.com/LadInTheLab/repairpbc/releases)
2. Download the executable for your operating system:
   - **Windows**: `repairpbc-windows.exe`
   - **macOS**: `repairpbc-macos`
   - **Linux**: `repairpbc-linux`
3. Make it executable (Linux/macOS): `chmod +x repairpbc-macos`
4. Run: `./repairpbc-macos -t trajectory.xtc -p topology.gro -o output.xtc`

**Note**: Executables are approximately 80MB due to included dependencies.

### Option 2: Python Package (Recommended for Developers)
For users comfortable with Python or who want to modify the code:

1. Clone or download this repository
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run: `python repair_pbc.py -t trajectory.xtc -p topology.gro -o output.xtc`

## Usage

### Basic Usage
```bash
python repair_pbc.py --trajectory your_trajectory.xtc --topology your_topology.gro --output repaired_trajectory.xtc
```

### Command Line Options
- `--trajectory` or `-t`: Input trajectory file (.xtc format)
- `--topology` or `-p`: Input topology file (.gro format)
- `--output` or `-o`: Output file path (extension is set automatically for custom formats)
- `--output-mode`: Output mode (`full`, `protein-only`, or `custom-format`)
- `--output-format`: Output format for custom-format mode (`xtc`, `pdb`, `xyz`, `npy`)
- `--no-verbose`: Disable verbose output
- `--help` or `-h`: Show help message

### Examples

Output only the protein atoms:
```bash
python repair_pbc.py -t simulation.xtc -p structure.gro -o repaired.xtc --output-mode protein-only
```

Output corrected coordinates as a multi-frame PDB:
```bash
python repair_pbc.py -t simulation.xtc -p structure.gro -o repaired.pdb --output-mode custom-format --output-format pdb
```

Output corrected coordinates as a NumPy array:
```bash
python repair_pbc.py -t simulation.xtc -p structure.gro -o repaired.npy --output-mode custom-format --output-format npy
```

## Output Modes

- `full` (default):  
  Overwrites only the protein atom positions in the full trajectory. All other atoms (water, ions, ligands) are left unchanged. This may result in physically invalid structures (overlaps, broken waters, etc.). Use this mode if you want to keep the original atom count and topology, but be aware of the risks.

- `protein-only`:  
  Outputs a new `.xtc` and `.gro` file containing only the protein atoms. All other atoms are removed. This is the safest mode for analysis or visualization of just the protein. The resulting files cannot be used for full-system MD restarts.

- `custom-format`:  
  Outputs only the corrected protein coordinates in a user-specified format. Use the `--output-format` option to choose the format: `xtc`, `pdb`, `xyz`, or `npy`. This is useful for downstream analysis or visualization in other tools.

## Warnings

- The `full` mode may create physically invalid structures. Use with caution.
- The `protein-only` and `custom-format` modes are recommended for most analysis and visualization tasks.
- Always check your output files before using them in downstream applications.

## How it works

1. Validates input files to ensure a single, contiguous protein chain.
2. Loads your .gro topology file and .xtc trajectory file.
3. Identifies protein atoms (amino acids only).
4. Tests for periodic boundary condition "jumps" by ensuring that no box-sized translations of the atoms *reduce* inter-atomic distance - if they do, they are applied and propagated to downstream atoms. 


## Validation

The tool automatically validates your input files and will stop with clear error messages if:

- Non-contiguous chains detected: Gaps in residue numbering.
- Multiple chains detected: Large jumps in residue numbers (>100 residues).
- No protein atoms found: Only water/solvent in the system.
- File format errors: Invalid .gro or .xtc files.
- Mismatched data: Trajectory and topology don't match.

## Requirements

- Python 3.7 or higher
- GROMACS trajectory files (.xtc)
- GROMACS topology files (.gro)
- Sufficient disk space for output file

## Notes

- The tool only processes protein atoms (amino acids), ignoring water and other molecules.
- Processing time depends on trajectory size and number of frames.
- This tool requires cubic or orthogonal simulation boxes. Triclinic boxes (with non-90° angles) are not supported.
- Original files are never modified - always creates new output files.

## Troubleshooting

**File not found errors:**
- Ensure your .xtc and .gro files exist and are in the correct directory.
- Check file permissions.

**Non-contiguous protein chain errors:**
- Your system contains multiple protein chains or gaps in residue numbering.
- This tool only works with single, contiguous polypeptide chains.
- Consider using a different tool for multi-chain systems or pre-processing your files.

**Multiple chains detected errors:**
- Large jumps in residue numbers indicate multiple protein molecules.
- Ensure your .gro file contains only one protein chain.
- Check if your system has multiple chains (A, B, C, etc.) that need to be processed separately.

**Memory errors:**
- Large trajectories may require significant RAM.

**Format errors:**
- Ensure files are valid GROMACS format.
- Check that .gro and .xtc files are compatible (same system).

## License

This tool is provided as-is for academic and research use. 