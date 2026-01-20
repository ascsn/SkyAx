# Automated Blocking Selection for ¹²¹Sn

## Quick Start

```bash
# Test with 2 states (fast)
python3 quick_block.py for005 --test

# Full scan with 6 states, 8 CPUs
python3 quick_block.py for005 --cpus 8 --max-states 6

# Custom parallelization: 4 processes × 2 threads = 8 total
python3 auto_block.py --input for005 \
    --exec ../../skyax.gfortran.omp --processes 4 --threads 2 --max-states 6
```

## How It Works

1. Runs ¹²¹Sn with `istblo=(0,0)` to get single-particle spectrum
2. Finds neutron states near Fermi surface (0.1 < v² < 0.99)
3. Tests each blocking candidate in parallel
4. Returns lowest energy configuration

## Results

For ¹²¹Sn (50p + 71n, odd-N):

```
Neutron candidates:
  State 477: K= 1, jz= 5.5, v²=0.317, E= -7.37 MeV
  State 469: K=-1, jz= 4.5, v²=0.249, E= -7.34 MeV
  State 398: K=-1, jz= 1.5, v²=0.833, E= -8.12 MeV  ← BEST
  State 457: K= 1, jz= 3.5, v²=0.166, E= -7.32 MeV

Recommended: istblo=(0,398)
Energy: -1024.399 MeV
```

Best configuration saved to `autoblock/for005_best`

## Parallelization

Total CPUs = `processes` × `threads`

- **processes**: Number of states tested simultaneously (Python level)
- **threads**: OpenMP threads per calculation

Examples:
- `--processes 4 --threads 2`: Test 4 states at once, 2 threads each (8 CPUs)
- `--processes 2 --threads 4`: Test 2 states at once, 4 threads each (8 CPUs)
- `--processes 1 --threads 8`: Serial state testing, 8 threads per state (8 CPUs)

Use more processes when testing many candidates, more threads for larger nuclei.

## Output Files

```
autoblock/
├── initial_out_details.sk2    # Single-particle spectrum (no blocking)
├── block_398_out_details.sk2  # Results for each tested state
├── for005_best                # Input with best istblo setting
```

## Options

```
--cpus N          Total CPUs to use (auto-configures processes/threads)
--max-states N    Number of candidates to test (default: 6)
--test            Quick test mode (2 states only)
--workdir DIR     Working directory (default: autoblock)
```

## Notes

- Conjugate pairs filtered (only one per pair tested)
- Prefers positive j_z when selecting from conjugates
- State with v² ≈ 0.5 is at Fermi surface
- State with v² ≈ 0.8-0.9 is occupied (good for blocking)
