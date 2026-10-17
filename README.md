# latticeboltzmann-python

> **Calendar reconstruction:** The historical records in [`docs/reconstructed-calendar-2026.csv`](docs/reconstructed-calendar-2026.csv) were generated on 2026-10-18. Their historical author dates are synthetic; they do not mean this project, its improvements, or work on them existed on those dates. The actual software improvements were completed on 2026-10-18. Original upstream authorship and history are preserved. Calendar coverage is not a measure of engineering quality or historical working time.

Development edition based on [Philip Mocz's original project](https://github.com/pmocz/latticeboltzmann-python).
Original authorship, Git history, and the GPL-3.0 license are retained.

## Local improvement and verification

The D2Q9 equilibrium calculation reuses its population buffer and computes squared
fluid speed once per step instead of once per lattice direction. Three numerical
tests check density and momentum moments, the isotropic pressure tensor at rest,
and full overwriting of a reused buffer.

```bash
python -m pip install -r requirements.txt
python -m unittest discover -v
```

Verified on 2026-10-18 with Python 3.13.5 and NumPy 2.1.3. A 64x32, 25-step
simulation produced exactly the same population array as upstream. At 400x100,
a local equilibrium-kernel benchmark measured 20.87 ms before and 13.79 ms after
(median of seven repeats of twenty calls). This is a kernel measurement;
whole-simulation speedups and timings depend on the workload and machine.

---

Lattice Boltzmann fluid simulation

## Create Your Own Lattice Boltzmann Simulation (With Python)

### Philip Mocz (2020) Princeton University, [@PMocz](https://twitter.com/PMocz)

### [📝 Read the Algorithm Write-up on Medium](https://medium.com/swlh/create-your-own-lattice-boltzmann-simulation-with-python-8759e8b53b1c)

Simulate a flow past cylinder with the Lattice Boltzmann method


```
python latticeboltzmann.py
```

![Simulation](./latticeboltzmann.png)
