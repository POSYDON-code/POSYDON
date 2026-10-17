# Example grid ini files

Ini files for `posydon-grid`, in the new format, filled in for the v2.0
physics at Z = 0.0142. Copy one next to your grid's `grid.csv`, fill in the
`<...>` values (the SLURM account and partition) and run

```bash
posydon-grid explain HMS-HMS.ini   # where every inlist value comes from
posydon-grid setup HMS-HMS.ini     # write the run directory
```

| File | Grid |
|---|---|
| [HMS-HMS.ini](HMS-HMS.ini) | two hydrogen-rich stars from ZAMS |
| [CO-HMS_RLO.ini](CO-HMS_RLO.ini) | a compact object and a hydrogen-rich star from ZAMS |
| [CO-HeMS.ini](CO-HeMS.ini) | a compact object and a helium star built in formation steps |
| [HeMS-HMS.ini](HeMS-HMS.ini) | a helium star (formation steps) and a hydrogen-rich star from ZAMS |
| [single_HMS.ini](single_HMS.ini) | single hydrogen-rich stars from ZAMS |
| [single_HeMS.ini](single_HeMS.ini) | single helium stars, built in formation steps |
| [HMS-HMS_rerun_and_my_physics.ini](HMS-HMS_rerun_and_my_physics.ini) | a rerun of v2.0, with your own inlists and Fortran source on top |
| [HMS-HMS_dynamic.ini](HMS-HMS_dynamic.ini) | a dynamic grid: psy-cris picks the points, as one MPI job |
| [HMS-HMS_laptop.ini](HMS-HMS_laptop.ini) | a few points on a laptop or workstation, without SLURM |

They read the inlist repository from `~/POSYDON-MESA-INLISTS`, a clone of
[POSYDON-MESA-INLISTS](https://github.com/POSYDON-code/POSYDON-MESA-INLISTS)
on its `v2_sources` branch:

```bash
git clone --branch v2_sources https://github.com/POSYDON-code/POSYDON-MESA-INLISTS.git ~/POSYDON-MESA-INLISTS
```

The examples show the options that matter for each grid. For every option
of a grid type, with a note on each, write a starter file:

```bash
posydon-grid init CO-HeMS full.ini
```

The ini format is described in the POSYDON documentation, "The ini file"
(`docs/_source/components-overview/mesa_grids/inifile.rst`). These files are
read by POSYDON's unit tests (`test_setup_examples.py`), so they stay in step
with the format.
