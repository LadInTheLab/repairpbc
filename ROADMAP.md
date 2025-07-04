# RepairPBC Development Roadmap

## Current Version (v1.0.0)
- Basic PBC repair for single protein chains
- Support for cubic/orthogonal simulation boxes
- Multiple output formats (xtc, pdb, xyz, npy)
- Comprehensive validation and error handling
- Cross-platform executables

## Planned Features

### v1.1.0 (Short-term)
- Enhanced validation
  - Warning system for potential issues

### v1.2.0 (Medium-term)
- Triclinic box support
  - Coordinate transformations for non-90° angles
  - Support for skewed simulation boxes
  - Validation of triclinic box parameters

### v2.0.0 (Long-term)
- Multi-chain support
  - Processing multiple protein chains
  - Chain-specific PBC repair
  - Inter-chain distance validation

## Future Considerations

### Research Applications
- Integration with analysis pipelines
- Batch processing of multiple trajectories

### User Experience
- Graphical user interface (GUI)
- Integration with molecular visualization software

## Contributing

I welcome contributions. If you'd like to work on any of these features:

1. Check existing issues for similar proposals
2. Create a new issue describing your proposed feature
3. Fork the repository and create a feature branch
4. Submit a pull request with your implementation

## Notes

- Priority is given to features that improve the core PBC repair functionality
- Performance improvements are always welcome
- New features should maintain backward compatibility
- All changes must include appropriate tests and documentation

---

This roadmap is a living document and will be updated as the project evolves. 