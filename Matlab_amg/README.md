# AMG-Preconditioned PCG for Pore-Network Models

AMG+ichol preconditioned conjugate gradient solver for SPD systems.

## File Structure

```
src/
    ├── test_amg.m                  # Main test driver
    ├── compile_mex.m               # Compiles all MEX files
    │
    ├── amg_setup.m                 # Build AMG hierarchy + fine-level ichol
    ├── amg_precondition.m          # Apply AMG+ichol preconditioner
    ├── amg_backslash.m             # Descending half-cycle
    ├── amg_fwdslash.m              # Ascending half-cycle
    │
    ├── rs_strength.m               # Wrapper → rs_strength_mex
    ├── rs_strength_mex.c           # Strength of connection (C/MEX)
    ├── rs_coarsening.m             # Wrapper → rs_coarsening_mex
    ├── rs_coarsening_mex.c         # Ruge-Stuben C/F splitting (C/MEX)
    ├── rs_interpolation.m          # Wrapper → rs_interpolation_mex
    ├── rs_interpolation_mex.c      # Direct interpolation (C/MEX)
    ├── cf_gauss_seidel.m           # Wrapper → cf_gauss_seidel_mex
    ├── cf_gauss_seidel_mex.c       # CF-ordered Gauss-Seidel smoother (C/MEX)
    │
    ├── read_nwk.m                  # Read raw network files and build system matrix
    ├── caseReaderMJ.m              # Network file reader
    ├── ConnectivityMx.m            # Build vertex-edge incidence matrix
    ├── makeDiffusiveConductance.m  # Compute edge conductances
    │
    └── four_large_Cartesian_networks/  # Test data (place here)
```

## Start

1. Place `four_large_Cartesian_networks/` inside `src/`.

2. **Compile MEX files** :
   ```matlab
   >> cd src
   >> compile_mex
   ```

3. **Run tests**:
   ```matlab
   >> test_amg                                              % default: 100x100x10 + 100x100x100
   >> test_amg({'100x100x10'})                              % single size
   >> test_amg({'100x100x10', '100x100x100', '100x100x1000'})  % multiple sizes
   ```

   Available test cases: `100x100x10`, `100x100x100`, `100x100x1000`, `100x1000x1000`.

## Requirements

- MATLAB 
- A C compiler
