# RepairPBC v[VERSION]

## What's New
- [List new features/fixes]

## Downloads

### For End Users (No Python Required)
Download the executable for your operating system:

- **Windows**: `repairpbc-windows.exe` - Run directly on Windows
- **macOS**: `repairpbc-macos` - Run on macOS (Intel/Apple Silicon)
- **Linux**: `repairpbc-linux` - Run on Linux systems

**Quick Start *(For those with some experience):***
1. Download the file for your OS
2. Make executable (Linux/macOS): `chmod +x repairpbc-[os]`
3. Run: `./repairpbc-[os] -t trajectory.xtc -p topology.gro -o output.xtc`

**Quick Start *(If you're new to command line tools):***
1. Download the file for your OS
2. *Optional: put the file in a folder that will be your "working directory"*
3. Open a command line window (cmd on Windows, Terminal on Mac)
3. Navigate to the location where the downloaded file is located using `cd`, as in `cd Downloads` if the file is in 
Downloads
2. Make executable (Linux/macOS): type`chmod +x repairpbc-[os]` and hit enter - on Mac, this would look like `chmod +x repairpbc-macos`
3. Run: `./repairpbc-[os] -t trajectory.xtc -p topology.gro -o output.xtc`

### For Developers
- Clone the repository and run `pip install -r requirements.txt`
- See the README for full installation instructions

## System Requirements
- Windows 10+, macOS 10.14+, or Linux
- 4GB+ RAM recommended for large trajectories
- No Python installation required for executables

## Notes
- Executables are ~80MB due to included dependencies
- Works with GROMACS .xtc and .gro files
- Only processes single, contiguous protein chains 