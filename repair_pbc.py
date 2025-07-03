#!/usr/bin/env python3
"""
RepairPBC - GROMACS Trajectory PBC Repair Tool

A command-line tool to repair periodic boundary conditions in GROMACS trajectory files.
"""

import os
import sys
import time
import argparse
import numpy as np
import mdtraj as mdt
import pandas as pd
from pathlib import Path

# Amino acid definitions
AMINO_ACIDS = [
    ("Alanine", "Ala", "A"), ("Arginine", "Arg", "R"), ("Asparagine", "Asn", "N"),
    ("Aspartic Acid", "Asp", "D"), ("Cysteine", "Cys", "C"), ("Glutamine", "Gln", "Q"),
    ("Glutamic Acid", "Glu", "E"), ("Glycine", "Gly", "G"), ("Histidine", "His", "H"),
    ("Isoleucine", "Ile", "I"), ("Leucine", "Leu", "L"), ("Lysine", "Lys", "K"), 
    ("Methionine", "Met", "M"), ("Phenylalanine", "Phe", "F"), ("Proline", "Pro", "P"), 
    ("Serine", "Ser", "S"), ("Threonine", "Thr", "T"), ("Tryptophan", "Trp", "W"), 
    ("Tyrosine", "Tyr", "Y"), ("Valine", "Val", "V")
]

# GROMACS .gro file field definitions
GRO_FIELDS = [
    ["Residue_Num", int, 0, 5], ["Residue_Name", str, 5, 5], ["Atom_Name", str, 10, 5],
    ["Atom_Num", int, 15, 5], ["Position_X", float, 20, 8], ["Position_Y", float, 28, 8],
    ["Position_Z", float, 36, 8], ["Velocity_X", float, 44, 8], ["Velocity_Y", float, 52, 8],
    ["Velocity_Z", float, 60, 8]
]


class RepairPBCError(Exception):
    """Custom exception for RepairPBC errors."""
    pass


def validate_file_path(file_path, file_type, must_exist=True):
    """Validate file path and extension."""
    if not file_path:
        raise RepairPBCError(f"{file_type} file path cannot be empty")
    
    path = Path(file_path)
    
    if must_exist and not path.exists():
        raise RepairPBCError(f"{file_type} file not found: {file_path}")
    
    if not path.suffix.lower() == file_type.lower():
        raise RepairPBCError(f"{file_type} file must have .{file_type} extension: {file_path}")
    
    return path


def validate_output_path(output_path, output_mode='full', output_format=None):
    """Validate output file path and ensure directory exists."""
    if not output_path:
        raise RepairPBCError("Output file path cannot be empty")
    
    path = Path(output_path)
    
    # Ensure output directory exists
    path.parent.mkdir(parents=True, exist_ok=True)
    
    # For custom-format mode, check extension compatibility
    if output_mode == 'custom-format':
        if output_format is None:
            raise RepairPBCError("You must specify --output-format for custom-format mode")
        
        expected_ext = f'.{output_format}'
        current_ext = path.suffix.lower()
        
        if current_ext == '':
            # No extension provided, that's fine - we'll add it
            return path
        elif current_ext == expected_ext:
            # Correct extension provided, use as-is
            return path
        else:
            # Wrong extension provided
            raise RepairPBCError(
                f"Output file has wrong extension. Expected {expected_ext} for {output_format} format, "
                f"but got {current_ext}. Use -o filename{expected_ext} or just -o filename"
            )
    
    # For other modes, expect .xtc extension
    if not path.suffix.lower() == '.xtc':
        raise RepairPBCError(f"Output file must have .xtc extension for {output_mode} mode: {output_path}")
    
    return path


def validate_protein_chain_contiguity(topology_df, verbose=True):
    """
    Validate that protein atoms form a single, contiguous chain.
    
    This function checks for:
    1. Multiple protein chains (different molecules)
    2. Non-contiguous protein segments (gaps in residue numbering)
    3. Chain boundaries (where residue numbers jump or reset)
    
    Args:
        topology_df: DataFrame containing protein topology information
        verbose: Enable verbose output
        
    Returns:
        bool: True if valid single contiguous chain, False otherwise
        
    Raises:
        RepairPBCError: If validation fails with specific error message
    """
    if len(topology_df) == 0:
        raise RepairPBCError("No protein atoms found in topology")
    
    # Get residue numbers and sort them
    residue_nums = topology_df["TOPOLOGY_Residue_Num"].values
    unique_residues = np.unique(residue_nums)
    
    if verbose:
        print(f"Found {len(unique_residues)} unique residues with {len(topology_df)} total protein atoms")
    
    # Check for gaps in residue numbering
    sorted_residues = np.sort(unique_residues)
    residue_diffs = np.diff(sorted_residues)
    
    # Look for gaps larger than 1 (indicating non-contiguous segments)
    gaps = residue_diffs > 1
    gap_positions = np.where(gaps)[0]
    
    if len(gap_positions) > 0:
        gap_info = []
        for pos in gap_positions:
            before_gap = sorted_residues[pos]
            after_gap = sorted_residues[pos + 1]
            gap_size = after_gap - before_gap - 1
            gap_info.append(f"residue {before_gap} to {after_gap} (gap of {gap_size} residues)")
        
        raise RepairPBCError(
            f"Non-contiguous protein chain detected! Gaps found between: {', '.join(gap_info)}. "
            f"This tool only works with single, contiguous protein chains. "
            f"Please ensure your protein structure contains only one continuous polypeptide chain."
        )
    
    # Check for potential multiple chains by looking for large jumps in residue numbers
    # (Some systems might have residue numbers that wrap around or are very large)
    if len(sorted_residues) > 1:
        max_jump = np.max(residue_diffs)
        if max_jump > 100:  # Arbitrary threshold - adjust as needed
            raise RepairPBCError(
                f"Large jump in residue numbering detected (max jump: {max_jump}). "
                f"This may indicate multiple protein chains or non-standard residue numbering. "
                f"This tool requires a single, contiguous protein chain."
            )
    
    # Check if atoms within each residue are properly ordered
    # (This is a basic check - more sophisticated bond-based validation could be added)
    for residue_num in unique_residues:
        residue_atoms = topology_df[topology_df["TOPOLOGY_Residue_Num"] == residue_num]
        atom_nums = residue_atoms["TOPOLOGY_Atom_Num"].values
        
        # Check if atom numbers are sequential within the residue
        if len(atom_nums) > 1:
            atom_diffs = np.diff(np.sort(atom_nums))
            if not np.all(atom_diffs == 1):
                if verbose:
                    print(f"Warning: Non-sequential atom numbering in residue {residue_num}")
    
    if verbose:
        print(f"Protein chain validation passed: {len(unique_residues)} contiguous residues")
    
    return True


def validate_trajectory_compatibility(trajectory_array, topology_df, verbose=True):
    """
    Validate that trajectory data is compatible with the topology.
    
    Args:
        trajectory_array: 3D numpy array of trajectory coordinates
        topology_df: DataFrame containing protein topology information
        verbose: Enable verbose output
        
    Returns:
        bool: True if compatible, False otherwise
    """
    if len(trajectory_array.shape) != 3:
        raise RepairPBCError(f"Trajectory must be 3D array, got shape: {trajectory_array.shape}")
    
    if trajectory_array.shape[2] != 3:
        raise RepairPBCError(f"Trajectory must contain 3D coordinates, got shape: {trajectory_array.shape}")
    
    # Check that number of atoms in trajectory matches topology
    n_atoms_traj = trajectory_array.shape[0]
    n_atoms_top = len(topology_df)
    
    if n_atoms_traj != n_atoms_top:
        raise RepairPBCError(
            f"Atom count mismatch: trajectory has {n_atoms_traj} atoms, "
            f"topology has {n_atoms_top} protein atoms"
        )
    
    if verbose:
        print(f"✓ Trajectory compatibility validated: {n_atoms_traj} atoms, {trajectory_array.shape[1]} frames")
    
    return True


def gro2df(gro_path, verbose=True):
    """Load and parse GROMACS .gro topology file."""
    try:
        if verbose:
            print(f"Loading topology file: {gro_path}")
        
        with open(gro_path, 'r') as fid:
            lines = fid.readlines()
        
        if len(lines) < 3:
            raise RepairPBCError(f"Invalid .gro file: insufficient lines in {gro_path}")
        
        title, num_atoms_line = lines[0], lines[1]
        
        try:
            num_atoms = int(num_atoms_line.strip())
        except ValueError:
            raise RepairPBCError(f"Invalid .gro file: cannot parse number of atoms in {gro_path}")
        
        # Get the lines containing atom data
        data_lines = lines[2:len(lines) - 1]
        
        if len(data_lines) != num_atoms:
            raise RepairPBCError(f"Atom count mismatch in {gro_path}: expected {num_atoms}, found {len(data_lines)}")
        
        # Parse atom data
        data_line_lists = []
        for i, line in enumerate(data_lines):
            if len(line) < 68:  # Minimum length for .gro format
                raise RepairPBCError(f"Invalid atom line {i+3} in {gro_path}: line too short")
            
            try:
                parsed_line = [line[field[2]:field[2] + field[3]].strip() for field in GRO_FIELDS]
                data_line_lists.append(parsed_line)
            except IndexError:
                raise RepairPBCError(f"Invalid atom line {i+3} in {gro_path}: incorrect format")
        
        # Build dataframe
        dataframe = pd.DataFrame(data_line_lists, columns=[f"TOPOLOGY_{field[0]}" for field in GRO_FIELDS])
        
        # Convert numeric columns
        numeric_cols = [f"TOPOLOGY_{field[0]}" for field in GRO_FIELDS if field[1] != str]
        dataframe[numeric_cols] = dataframe[numeric_cols].apply(pd.to_numeric, errors='coerce')
        
        # Check for conversion errors
        if dataframe[numeric_cols].isnull().any().any():
            raise RepairPBCError(f"Invalid numeric data in {gro_path}")
        
        # Extract protein atoms only
        protein_atoms = [aa[1].upper() for aa in AMINO_ACIDS]
        selected_topology_df = dataframe[dataframe["TOPOLOGY_Residue_Name"].isin(protein_atoms)]
        
        if len(selected_topology_df) == 0:
            raise RepairPBCError(f"No protein atoms found in {gro_path}")
        
        # Parse box dimensions
        box_line = lines[-1].strip()
        try:
            box_vec_list = box_line.split()
            if len(box_vec_list) < 3:
                raise ValueError("Insufficient box dimensions")
            box_dims = np.array([float(box_vec_list[0]), float(box_vec_list[1]), float(box_vec_list[2])])
        except (ValueError, IndexError):
            raise RepairPBCError(f"Invalid box dimensions in {gro_path}")
        
        if verbose:
            print(f"Topology loaded: {len(selected_topology_df)} protein atoms, box dimensions: {box_dims}")
        
        return selected_topology_df, box_dims, dataframe
        
    except FileNotFoundError:
        raise RepairPBCError(f"Topology file not found: {gro_path}")
    except Exception as e:
        raise RepairPBCError(f"Error loading topology file {gro_path}: {str(e)}")


def xtc2df(xtc_path, topology_path, verbose=True):
    """Load and parse GROMACS .xtc trajectory file."""
    try:
        if verbose:
            print(f"Loading trajectory file: {xtc_path}")
        
        # Load topology info
        structure_topology_df, box_dims, _ = gro2df(topology_path, verbose=False)
        
        # Load trajectory using MDTraj
        try:
            trajectory_object = mdt.load_xtc(xtc_path, topology_path)
        except Exception as e:
            raise RepairPBCError(f"Error loading trajectory with MDTraj: {str(e)}")
        
        if trajectory_object.n_frames == 0:
            raise RepairPBCError(f"Trajectory file {xtc_path} contains no frames")
        
        # Extract protein atom coordinates
        selected_topology_indices = structure_topology_df.index
        selected_trajectories = trajectory_object.xyz[:, selected_topology_indices, :]
        selected_trajectories = np.swapaxes(selected_trajectories, 0, 1)
        
        if verbose:
            print(f"Trajectory loaded: {trajectory_object.n_frames} frames, {selected_trajectories.shape[1]} protein atoms")
        
        return selected_trajectories, structure_topology_df, box_dims
        
    except Exception as e:
        raise RepairPBCError(f"Error loading trajectory file {xtc_path}: {str(e)}")


def repair_pbc(trajectory_file_path, topology_file_path, verbose=True):
    """
    Repair periodic boundary conditions in GROMACS trajectory.
    
    This function requires a single, contiguous protein chain. It will validate
    the input files to ensure they meet this requirement before processing.
    
    Args:
        trajectory_file_path: Path to .xtc trajectory file
        topology_file_path: Path to .gro topology file
        verbose: Enable verbose output
        
    Returns:
        corrected_trajectory: numpy array of corrected coordinates
        
    Raises:
        RepairPBCError: If validation fails or processing encounters errors
    """
    try:
        if verbose:
            print("=" * 60)
            print("RepairPBC - GROMACS Trajectory PBC Repair Tool")
            print("=" * 60)
            print("IMPORTANT: This tool only works with single, contiguous protein chains.")
            print("Multiple chains or non-contiguous segments will cause incorrect results.")
            print("-" * 60)
        
        # Load files
        structure_topology_df, box_dimensions, _ = gro2df(topology_file_path, verbose=verbose)
        
        # Validate protein chain contiguity BEFORE loading trajectory
        validate_protein_chain_contiguity(structure_topology_df, verbose=verbose)
        
        # Load trajectory
        selected_trajectories, _, _ = xtc2df(trajectory_file_path, topology_file_path, verbose=verbose)
        selected_trajectories_array = np.array(selected_trajectories)
        
        # Validate trajectory compatibility
        validate_trajectory_compatibility(selected_trajectories_array, structure_topology_df, verbose=verbose)
        
        # Get box parameters
        offsets = box_dimensions
        box_center = offsets / 2
        allowed_multipliers = np.arange(0, 5) * 0.5
        
        start_time = time.perf_counter()
        if verbose:
            print(f"Starting PBC correction...")
        
        # Calculate allowed corrections
        allowed_corrections_X = allowed_multipliers * offsets[0]
        allowed_corrections_Y = allowed_multipliers * offsets[1]
        allowed_corrections_Z = allowed_multipliers * offsets[2]
        
        # Calculate displacements
        trajectory_upshifted = np.roll(selected_trajectories_array, shift=-1, axis=0)
        fwd_displacements = trajectory_upshifted - selected_trajectories_array
        
        # Find correction signs
        fwd_displacements_signs = np.sign(fwd_displacements)
        fwd_correction_signs = fwd_displacements_signs * -1
        
        # Broadcast corrections
        allowed_corrections_X_broadcast = np.tile(
            allowed_corrections_X, (fwd_displacements.shape[0], fwd_displacements.shape[1], 1)
        )
        allowed_corrections_Y_broadcast = np.tile(
            allowed_corrections_Y, (fwd_displacements.shape[0], fwd_displacements.shape[1], 1)
        )
        allowed_corrections_Z_broadcast = np.tile(
            allowed_corrections_Z, (fwd_displacements.shape[0], fwd_displacements.shape[1], 1)
        )
        
        # Apply signs to corrections
        allowed_corrections_X_signed = allowed_corrections_X_broadcast * fwd_correction_signs[:, :, 0][:, :, np.newaxis]
        allowed_corrections_Y_signed = allowed_corrections_Y_broadcast * fwd_correction_signs[:, :, 1][:, :, np.newaxis]
        allowed_corrections_Z_signed = allowed_corrections_Z_broadcast * fwd_correction_signs[:, :, 2][:, :, np.newaxis]
        
        # Apply corrections to displacements
        corrected_disp_X = fwd_displacements[:, :, 0][:, :, np.newaxis] + allowed_corrections_X_signed
        corrected_disp_Y = fwd_displacements[:, :, 1][:, :, np.newaxis] + allowed_corrections_Y_signed
        corrected_disp_Z = fwd_displacements[:, :, 2][:, :, np.newaxis] + allowed_corrections_Z_signed
        
        # Find minimum displacements
        idx_min_displacement_X = np.argmin(np.abs(corrected_disp_X), axis=2)
        idx_min_displacement_Y = np.argmin(np.abs(corrected_disp_Y), axis=2)
        idx_min_displacement_Z = np.argmin(np.abs(corrected_disp_Z), axis=2)
        
        # Extract best corrections
        best_X_corrections = allowed_corrections_X_signed[
            np.arange(allowed_corrections_X_signed.shape[0])[:, None],
            np.arange(allowed_corrections_X_signed.shape[1]), idx_min_displacement_X
        ]
        best_Y_corrections = allowed_corrections_Y_signed[
            np.arange(allowed_corrections_Y_signed.shape[0])[:, None],
            np.arange(allowed_corrections_Y_signed.shape[1]), idx_min_displacement_Y
        ]
        best_Z_corrections = allowed_corrections_Z_signed[
            np.arange(allowed_corrections_Z_signed.shape[0])[:, None],
            np.arange(allowed_corrections_Z_signed.shape[1]), idx_min_displacement_Z
        ]
        
        # Combine corrections
        best_corrections_XYZ = np.stack((best_X_corrections, best_Y_corrections, best_Z_corrections), axis=2)
        
        # Apply cumulative corrections
        best_corrections_shifted = np.roll(best_corrections_XYZ, shift=1, axis=0)
        cumulative_corrections = np.cumsum(best_corrections_shifted, axis=0)
        
        # Calculate corrected positions
        corrected_trajectory = selected_trajectories_array + cumulative_corrections
        
        # Center protein in box
        centroids = np.mean(corrected_trajectory, axis=0)
        centroids_shift = box_center - centroids
        centroids_shift_broadcast = np.broadcast_to(centroids_shift, corrected_trajectory.shape)
        corrected_trajectory = corrected_trajectory + centroids_shift_broadcast
        
        end_time = time.perf_counter()
        time_elapsed = np.round((end_time - start_time), decimals=2)
        total_points = (corrected_trajectory.shape[0] * corrected_trajectory.shape[1] * 
                       corrected_trajectory.shape[2]) + corrected_trajectory.shape[1]
        
        if verbose:
            print(f"PBC correction completed in {time_elapsed} seconds")
            print(f"Processed {total_points:,} inter-atomic distances")
        
        return corrected_trajectory
        
    except Exception as e:
        raise RepairPBCError(f"Error during PBC repair: {str(e)}")


def write_gro_file_protein_only(topology_path, output_gro_path, protein_indices, verbose=True):
    """
    Write a new .gro file containing only the protein atoms.
    """
    with open(topology_path, 'r') as fid:
        lines = fid.readlines()
    title = lines[0]
    num_atoms = int(lines[1].strip())
    data_lines = lines[2:len(lines) - 1]
    box_line = lines[-1]
    # Only keep lines for protein atoms
    protein_lines = [data_lines[i] for i in protein_indices]
    with open(output_gro_path, 'w') as fout:
        fout.write(title)
        fout.write(f"{len(protein_lines)}\n")
        for line in protein_lines:
            fout.write(line)
        fout.write(box_line)
    if verbose:
        print(f"Wrote new .gro file with only protein atoms: {output_gro_path}")


def save_trajectory(corrected_trajectory, output_path, topology_path, trajectory_path, output_mode='full', output_format=None, verbose=True):
    """Save corrected trajectory in the specified mode."""
    if verbose:
        if output_mode == 'full':
            print(f"Saving corrected trajectory to: {output_path}")
        elif output_mode == 'protein-only':
            print("Saving protein-only files...")
        elif output_mode == 'custom-format':
            if output_format == 'npy':
                if str(output_path).endswith('.npy'):
                    actual_out = str(output_path)
                else:
                    actual_out = str(output_path) + '.npy'
                print(f"Saving corrected coordinates as .npy: {actual_out}")
            elif output_format == 'xyz':
                if str(output_path).endswith('.xyz'):
                    actual_out = str(output_path)
                else:
                    actual_out = str(output_path) + '.xyz'
                print(f"Saving corrected coordinates as .xyz: {actual_out}")
            elif output_format == 'pdb':
                if str(output_path).endswith('.pdb'):
                    actual_out = str(output_path)
                else:
                    actual_out = str(output_path) + '.pdb'
                print(f"Saving corrected coordinates as .pdb: {actual_out}")
            else:
                print(f"Saving corrected trajectory to: {output_path}")
        else:
            print(f"Saving corrected trajectory to: {output_path}")
    if output_mode == 'full':
        print("Warning: Only protein atom positions are updated. Non-protein atoms may overlap with the protein or be physically invalid.")
        original_traj = mdt.load_xtc(str(trajectory_path), str(topology_path))
        protein_atoms = [aa[1].upper() for aa in AMINO_ACIDS]
        protein_indices = []
        for i, residue in enumerate(original_traj.topology.residues):
            if residue.name in protein_atoms:
                residue_atoms = original_traj.topology.select(f"resid {i}")
                protein_indices.extend(residue_atoms)
        full_xyz = original_traj.xyz.copy()
        corrected_xyz = np.swapaxes(corrected_trajectory, 0, 1)
        full_xyz[:, protein_indices, :] = corrected_xyz
        new_traj = mdt.Trajectory(xyz=full_xyz, topology=original_traj.topology)
        new_traj.save_xtc(str(output_path))
        if verbose:
            print(f"Trajectory saved in full mode: {output_path}")
    elif output_mode == 'protein-only':
        print("Outputting only the protein atoms. The resulting .xtc and .gro will not contain any other atoms.")
        # Load topology and get protein atom indices
        with open(topology_path, 'r') as fid:
            lines = fid.readlines()
        data_lines = lines[2:len(lines) - 1]
        protein_atoms = [aa[1].upper() for aa in AMINO_ACIDS]
        protein_indices = []
        for i, line in enumerate(data_lines):
            residue_name = line[5:10].strip().upper()
            if residue_name in protein_atoms:
                protein_indices.append(i)
        # Write new .gro file
        gro_out = str(output_path).replace('.xtc', '_protein.gro')
        write_gro_file_protein_only(topology_path, gro_out, protein_indices, verbose=verbose)
        # Save new .xtc with only protein atoms
        corrected_xyz = np.swapaxes(corrected_trajectory, 0, 1)
        # Create a minimal topology for the protein
        import tempfile
        temp_gro = tempfile.NamedTemporaryFile(delete=False, suffix='.gro')
        write_gro_file_protein_only(topology_path, temp_gro.name, protein_indices, verbose=False)
        protein_top = mdt.load(temp_gro.name)
        new_traj = mdt.Trajectory(xyz=corrected_xyz, topology=protein_top.topology)
        xtc_out = str(output_path).replace('.xtc', '_protein.xtc')
        new_traj.save_xtc(xtc_out)
        temp_gro.close()
        os.unlink(temp_gro.name)
        if verbose:
            print(f"Protein-only .xtc written to: {xtc_out}")
    elif output_mode == 'custom-format':
        if output_format is None:
            raise RepairPBCError("You must specify --output-format for custom-format mode.")
        corrected_xyz = np.swapaxes(corrected_trajectory, 0, 1)
        if output_format == 'npy':
            # Check if user already provided .npy extension
            if str(output_path).endswith('.npy'):
                npy_out = str(output_path)
            else:
                npy_out = str(output_path) + '.npy'
            np.save(npy_out, corrected_xyz)
            if verbose:
                print(f"Corrected coordinates saved as .npy: {npy_out}")
        elif output_format == 'xyz':
            # Save as multi-frame XYZ (simple format, no topology)
            n_frames, n_atoms, _ = corrected_xyz.shape
            if str(output_path).endswith('.xyz'):
                xyz_out = str(output_path)
            else:
                xyz_out = str(output_path) + '.xyz'
            with open(xyz_out, 'w') as fout:
                for frame in range(n_frames):
                    fout.write(f"{n_atoms}\nFrame {frame}\n")
                    for atom in range(n_atoms):
                        x, y, z = corrected_xyz[frame, atom]
                        fout.write(f"ATOM {x:.3f} {y:.3f} {z:.3f}\n")
            if verbose:
                print(f"Corrected coordinates saved as .xyz: {xyz_out}")
        elif output_format == 'pdb':
            # Save as multi-frame PDB (using MDTraj)
            import tempfile
            temp_gro = tempfile.NamedTemporaryFile(delete=False, suffix='.gro')
            with open(topology_path, 'r') as fid:
                lines = fid.readlines()
            data_lines = lines[2:len(lines) - 1]
            protein_atoms = [aa[1].upper() for aa in AMINO_ACIDS]
            protein_indices = []
            for i, line in enumerate(data_lines):
                residue_name = line[5:10].strip().upper()
                if residue_name in protein_atoms:
                    protein_indices.append(i)
            write_gro_file_protein_only(topology_path, temp_gro.name, protein_indices, verbose=False)
            protein_top = mdt.load(temp_gro.name)
            new_traj = mdt.Trajectory(xyz=corrected_xyz, topology=protein_top.topology)
            if str(output_path).endswith('.pdb'):
                pdb_out = str(output_path)
            else:
                pdb_out = str(output_path) + '.pdb'
            new_traj.save_pdb(pdb_out)
            temp_gro.close()
            os.unlink(temp_gro.name)
            if verbose:
                print(f"Corrected coordinates saved as .pdb: {pdb_out}")
        else:
            raise RepairPBCError(f"Unsupported output format: {output_format}")
    else:
        raise RepairPBCError(f"Unknown output mode: {output_mode}")


def main():
    """Main CLI function."""
    parser = argparse.ArgumentParser(
        description="Repair periodic boundary conditions in GROMACS trajectory files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python repair_pbc.py -t simulation.xtc -p structure.gro -o repaired.xtc
  python repair_pbc.py --trajectory traj.xtc --topology topol.gro --output fixed.xtc --no-verbose
        """
    )
    
    parser.add_argument(
        '-t', '--trajectory',
        required=True,
        help='Input trajectory file (.xtc format)'
    )
    
    parser.add_argument(
        '-p', '--topology', 
        required=True,
        help='Input topology file (.gro format)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output file path (extension will be set automatically based on output mode and format)'
    )
    
    parser.add_argument(
        '--no-verbose',
        action='store_true',
        help='Disable verbose output'
    )
    
    parser.add_argument(
        '--output-mode',
        choices=['full', 'protein-only', 'custom-format'],
        default='full',
        help='Output mode: full (default, overwrite only protein atoms in full trajectory), protein-only (output only protein atoms), custom-format (output in another format)'
    )
    
    parser.add_argument(
        '--output-format',
        choices=['xtc', 'pdb', 'xyz', 'npy'],
        default=None,
        help='Output format for custom-format mode (e.g., pdb, xyz, npy)'
    )
    
    args = parser.parse_args()
    
    # Set verbose flag
    verbose = not args.no_verbose
    
    try:
        # Validate input files
        trajectory_path = validate_file_path(args.trajectory, '.xtc')
        topology_path = validate_file_path(args.topology, '.gro')
        output_mode = args.output_mode
        output_format = args.output_format
        output_path = validate_output_path(args.output, output_mode, output_format)
        
        if verbose:
            print("=" * 60)
            print("RepairPBC - GROMACS Trajectory PBC Repair Tool")
            print("=" * 60)
            print(f"Input trajectory: {trajectory_path}")
            print(f"Input topology: {topology_path}")
            print(f"Output file: {output_path}")
            print("-" * 60)
        
        # Perform PBC repair
        corrected_trajectory = repair_pbc(trajectory_path, topology_path, verbose=verbose)
        
        # Save corrected trajectory
        save_trajectory(
            corrected_trajectory,
            output_path,
            topology_path,
            trajectory_path,
            output_mode=output_mode,
            output_format=output_format,
            verbose=verbose
        )
        
        if verbose:
            print("-" * 60)
            print("PBC repair completed successfully!")
            if output_mode == 'full':
                print(f"Output saved to: {output_path}")
            elif output_mode == 'protein-only':
                gro_out = str(output_path).replace('.xtc', '_protein.gro')
                xtc_out = str(output_path).replace('.xtc', '_protein.xtc')
                print(f"Protein-only files saved:")
                print(f"  Topology: {gro_out}")
                print(f"  Trajectory: {xtc_out}")
            elif output_mode == 'custom-format':
                if output_format == 'npy':
                    if str(output_path).endswith('.npy'):
                        actual_out = str(output_path)
                    else:
                        actual_out = str(output_path) + '.npy'
                    print(f"Output saved to: {actual_out}")
                elif output_format == 'xyz':
                    if str(output_path).endswith('.xyz'):
                        actual_out = str(output_path)
                    else:
                        actual_out = str(output_path) + '.xyz'
                    print(f"Output saved to: {actual_out}")
                elif output_format == 'pdb':
                    if str(output_path).endswith('.pdb'):
                        actual_out = str(output_path)
                    else:
                        actual_out = str(output_path) + '.pdb'
                    print(f"Output saved to: {actual_out}")
                else:
                    print(f"Output saved to: {output_path}")
            else:
                print(f"Output saved to: {output_path}")
            print("=" * 60)
        
    except RepairPBCError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"UNEXPECTED ERROR: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 