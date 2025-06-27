# Feature: SuPy Data Processing Enhancements

## Context
This branch focuses on improving data processing capabilities in SuPy, including better handling of input/output data, enhanced post-processing functions, and improved data validation and preprocessing workflows.

## GitHub Issues to Address
- **#408**: Full output from save_supy. Correct? (post-processing)
- **#417**: Output building ohm coefficients (OHM)
- **#412**: Update sample data with the benchmark one (pre-processing)
- **#406**: LAI and SMD not correct using sample data (pre-processing, user question)
- **#353**: Update STEBBS units to match SUEWS (STEBBS)

## Progress Tracking
- [ ] Fix save_supy output issues (#408)
  - [ ] Investigate full output correctness
  - [ ] Ensure all variables are properly saved
  - [ ] Fix any data loss or corruption
  - [ ] Improve output documentation
- [ ] Enhance OHM coefficient output (#417)
  - [ ] Add building-specific OHM coefficients to output
  - [ ] Create post-processing functions for OHM analysis
  - [ ] Document coefficient interpretation
- [ ] Improve data preprocessing (#406, #412)
  - [ ] Fix LAI and SMD preprocessing issues
  - [ ] Update sample data to use benchmark
  - [ ] Enhance data validation routines
  - [ ] Add data quality checks
- [ ] Unit standardisation (#353)
  - [ ] Ensure consistent units across SUEWS/STEBBS
  - [ ] Add unit conversion utilities
  - [ ] Document all unit conventions
- [ ] General data processing enhancements
  - [ ] Add more post-processing utilities
  - [ ] Improve data visualisation functions
  - [ ] Enhance gap-filling algorithms
  - [ ] Add data export formats

## Key Decisions
- Maintain backward compatibility for existing workflows
- Use standard data formats (NetCDF, HDF5) where appropriate
- Provide clear data provenance tracking
- Focus on user-friendly interfaces

## Implementation Notes
- Many issues stem from incomplete data handling
- Need better metadata preservation through processing chain
- Consider lazy loading for large datasets
- Implement proper data validation at each step

## Files to Modify
- `src/supy/_save.py` - Fix save_supy functionality
- `src/supy/_post.py` - Enhanced post-processing
- `src/supy/_load.py` - Improved data loading
- `src/supy/util/_io.py` - I/O utilities
- `src/supy/_check.py` - Data validation
- `src/supy/data_model/` - Data model enhancements

## New Features to Implement
1. **Enhanced Output System**
   - Complete variable output with metadata
   - Multiple output format support
   - Compressed output options
   - Streaming output for long runs

2. **Post-processing Toolkit**
   - OHM coefficient analysis tools
   - Energy balance closure functions
   - Footprint analysis utilities
   - Statistical summary functions

3. **Data Quality Tools**
   - Automatic quality flagging
   - Gap detection and filling
   - Outlier identification
   - Data consistency checks

4. **Visualisation Suite**
   - Time series plotting utilities
   - Spatial mapping functions
   - Energy balance diagrams
   - Comparative analysis plots

## Testing Approach
- Test with various input data formats
- Verify output completeness and accuracy
- Test with large datasets for performance
- Ensure unit consistency throughout