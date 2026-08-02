.. _memory-tuning:

Tuning Memory Usage
===================

A population run's memory use, and the related setting ``dump_rate``, control
the balance between speed and RAM. This page gives concrete guidance so you can
request the right amount of memory on an HPC facility and avoid extreme I/O.

How much memory do I need?
--------------------------

A population run needs a bare minimum of about **8 GB of memory per CPU**. The
DR2 grids have to be loaded into memory; binary grids are loaded at startup,
while single-star models are loaded on demand as they are needed.

However, 8 GB only lets you keep a small number of binaries in memory at once.
To stay within budget you would need a ``dump_rate`` below ~ 1000, which
forces many I/O operations and slows the simulation.

A better starting point is **9-10 GB per CPU**. This lets you keep more
binaries in memory (making the run faster) while still allowing a healthy
``dump_rate`` of a few thousand binaries.

What is ``dump_rate``?
-------------------------

``dump_rate`` is the number of binaries that are kept in memory before POSYDON
writes them out and frees that memory. It is a trade-off between RAM and I/O.

General guidance:

* use at least ``dump_rate`` **500** for populations of 100,000 binaries or
  more;
* a very low ``dump_rate`` on a large population forces many I/O calls (writing
  the binary data out and reading it back during the merging of output files)
  and can even introduce issues while reading, writing, and merging output;
* keep ``dump_rate`` relatively high, but not so high that it exceeds the number
  of binaries in the population (that would defeats the purpose, since you would
  only dump once at the very end).

Related
-------

* :ref:`run-hpc` for the walltime and array size that pair with a given
  ``dump_rate``.
