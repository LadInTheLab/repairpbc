#!/usr/bin/env python3
"""
Simple test script for RepairPBC tool.
This script tests the basic functionality without requiring actual GROMACS files.
"""

import sys
import os
from pathlib import Path

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from repair_pbc import validate_file_path, validate_output_path, RepairPBCError

def test_validation_functions():
    """Test the validation functions."""
    print("Testing validation functions...")
    
    # Test file path validation
    try:
        # This should fail - file doesn't exist
        validate_file_path("nonexistent.xtc", ".xtc")
        print("ERROR: Should have failed for nonexistent file")
        return False
    except RepairPBCError:
        print("Correctly caught nonexistent file error")
    
    try:
        # This should fail - wrong extension
        validate_file_path("test.txt", ".xtc")
        print("ERROR: Should have failed for wrong extension")
        return False
    except RepairPBCError:
        print("Correctly caught wrong extension error")
    
    # Test output path validation
    try:
        output_path = validate_output_path("test_output.xtc")
        print(f"Output path validation passed: {output_path}")
    except RepairPBCError as e:
        print(f"ERROR: Output path validation failed: {e}")
        return False
    
    print("All validation tests passed!")
    return True

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")
    
    try:
        import numpy as np
        print("NumPy imported successfully")
    except ImportError as e:
        print(f"ERROR: NumPy import failed: {e}")
        return False
    
    try:
        import pandas as pd
        print("Pandas imported successfully")
    except ImportError as e:
        print(f"ERROR: Pandas import failed: {e}")
        return False
    
    try:
        import mdtraj as mdt
        print("MDTraj imported successfully")
    except ImportError as e:
        print(f"ERROR: MDTraj import failed: {e}")
        return False
    
    print("All imports successful!")
    return True

def main():
    """Run all tests."""
    print("=" * 50)
    print("RepairPBC Tool Test Suite")
    print("=" * 50)
    
    tests = [
        ("Import Test", test_imports),
        ("Validation Test", test_validation_functions),
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        print(f"\nRunning {test_name}...")
        if test_func():
            passed += 1
            print(f"✓ {test_name} PASSED")
        else:
            print(f"✗ {test_name} FAILED")
    
    print("\n" + "=" * 50)
    print(f"Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("All tests passed! The tool is ready to use.")
        print("\nTo use the tool:")
        print("python repair_pbc.py -t your_trajectory.xtc -p your_topology.gro -o repaired.xtc")
    else:
        print("Some tests failed. Please check the errors above.")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 