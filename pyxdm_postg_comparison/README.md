# PyXDM vs postg

- **Integration grid:** postg uses a fixed mesh; PyXDM uses HORTON grids (`fine` by default).
- **Hirshfeld pro-atoms:** postg uses free-atom densities built into the program; PyXDM uses a pro-atom database at a chosen level of theory (here B3LYP/cc-pVTZ, `proatomdb/`).
- **Density cutoff:** postg sets b = 0 where the spin density is below 1e-10; PyXDM goes down to 1e-15.

## Pipeline

```bash
./build_postg.sh postg-build                                              # postg with postg.patch (dumps its mesh and weights)
python molecules.py geometries                                            # MMFF geometries
python run_orca.py geometries wfn --orca-dir /path/to/orca                # B3LYP/cc-pVTZ wavefunctions (molden)
python compare.py wfn results/comparison.csv --postg postg-build/postg    # postg; PyXDM on postg's mesh and weights; PyXDM defaults
python proatom_test.py wfn/*.molden.input --grid fine --output results/proatoms.csv \
    --atomicdata postg-build/atomicdata.f90                               # PyXDM with postg's pro-atoms
jupyter nbconvert --to notebook --execute --inplace analysis.ipynb        # tables and figures
```

The pro-atom database was made with `horton-atomdb.py input orca H,6-8,16 template.inp --max-cation 2 --max-anion 2`,
`bash run_orca.sh` and `horton-atomdb.py convert`.

## Note on ignored files

`work/`, which contains the mesh and grid files, is not uploaded due to its size (~1.4 GB) and must be regenerated using the pipeline above.
