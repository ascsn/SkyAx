#!/usr/bin/env python3
"""
Automated blocking state selection for SkyAx odd-nucleus calculations.

Physical principle:
- Run the target odd nucleus with NO blocking (istblo=0,0)
- Identify states near the Fermi surface (partial occupation v²~0.5)
- Test blocking each candidate state (in parallel at Python level)
- Select configuration with lowest energy

Features:
- Parses clean state format from bottom of out_details.sk2
- Filters conjugate pairs (prefers positive j_z)
- Python-level parallelization for testing multiple states simultaneously
- Proper file management to avoid conflicts
"""

import sys
import re
import subprocess
import shutil
from pathlib import Path
import argparse
import os
from multiprocessing import Pool, cpu_count
from functools import partial


class BlockingScanner:
    def __init__(self, input_file, executable, workdir='.', num_threads=1, num_processes=1):
        """
        Initialize the blocking scanner.
        
        Args:
            input_file: Path to input file (for TARGET odd nucleus)
            executable: Path to SkyAx executable
            workdir: Working directory for calculations
            num_threads: Number of OpenMP threads per process
            num_processes: Number of parallel Python processes
        """
        self.input_file = Path(input_file).resolve()
        self.executable = Path(executable).resolve()
        self.workdir = Path(workdir).resolve()
        self.num_threads = num_threads
        self.num_processes = num_processes
        
        self.workdir.mkdir(exist_ok=True, parents=True)
        
        if not self.input_file.exists():
            raise FileNotFoundError(f"Input file not found: {self.input_file}")
        if not self.executable.exists():
            raise FileNotFoundError(f"Executable not found: {self.executable}")
    
    def parse_input_nucleus(self, input_file):
        """Extract nprot and nneut from input file."""
        with open(input_file, 'r') as f:
            content = f.read()
        
        nprot_match = re.search(r'nprot\s*=\s*(\d+)', content)
        nneut_match = re.search(r'nneut\s*=\s*(\d+)', content)
        
        if not nprot_match or not nneut_match:
            raise ValueError("Could not find nprot/nneut in input file")
        
        return int(nprot_match.group(1)), int(nneut_match.group(1))
    
    def modify_input(self, input_file, output_file, istblo_p=None, istblo_n=None):
        """
        Modify input file to set blocking states.
        
        Args:
            input_file: Original input file path
            output_file: Path for modified input file
            istblo_p: Proton blocking state number (None to keep original)
            istblo_n: Neutron blocking state number (None to keep original)
        """
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Replace istblo if both values specified
        if istblo_p is not None and istblo_n is not None:
            pattern = r'istblo\s*=\s*\d+\s*,\s*\d+'
            replacement = f'istblo={istblo_p},{istblo_n}'
            
            if re.search(pattern, content):
                content = re.sub(pattern, replacement, content)
            else:
                # Add istblo to static namelist if not present
                content = re.sub(
                    r'(&static[^/]*)',
                    rf'\1  {replacement},\n',
                    content,
                    flags=re.DOTALL
                )
        
        with open(output_file, 'w') as f:
            f.write(content)
    
    def run_skyax(self, input_file, log_prefix):
        """
        Run SkyAx with given input file.
        
        Returns:
            float: Total energy in MeV (or None if failed)
        """
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(self.num_threads)
        
        log_file = self.workdir / f"{log_prefix}.log"
        
        try:
            with open(input_file, 'r') as fin, open(log_file, 'w') as fout:
                result = subprocess.run(
                    [str(self.executable)],
                    stdin=fin,
                    stdout=fout,
                    stderr=subprocess.STDOUT,
                    cwd=self.workdir,
                    env=env,
                    timeout=600  # 10 minute timeout
                )
            
            # Parse energy from output
            energy = None
            with open(log_file, 'r') as f:
                for line in f:
                    # Look for energy lines like "iter= 1000, energies: from rearr. vs. from densities= 0.0000 -1017.4357"
                    match = re.search(r'energies:\s+from\s+rearr.*=\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)', line)
                    if match:
                        energy = float(match.group(2))  # Use second value (from densities)
            
            # Archive output files with unique names
            for pattern in ['out_*']:
                for f in self.workdir.glob(pattern):
                    new_name = f.with_name(f"{log_prefix}_{f.name}")
                    if new_name.exists():
                        new_name.unlink()
                    shutil.move(str(f), str(new_name))
            
            return energy
            
        except subprocess.TimeoutExpired:
            print(f"    WARNING: {log_prefix} timed out")
            return None
        except Exception as e:
            print(f"    ERROR in {log_prefix}: {e}")
            return None
    
    def parse_state_details(self, details_file):
        """
        Parse out_details.sk2 file to extract single-particle state information.
        Uses the clean formatted section at the bottom of the file.
        
        Returns:
            list: List of state dictionaries with keys:
                  number, isospin, jz, K (tot), v2, energy, conjugate, etc.
        """
        states = []
        in_state_section = False
        
        with open(details_file, 'r') as f:
            for line in f:
                # Look for section headers
                if 'parities' in line and 'components' in line:
                    in_state_section = True
                    continue
                
                if not in_state_section:
                    continue
                
                # Parse state lines in format:
                # "  315 2 -0.5 -1 -1  2 1.00000 0.99994  0.01  0.53  -42.274  2.1720E-03   3.284      354"
                # Nr Iso j_z tot r z wtpair v^2 uv delta energy energy-variances s.p.radius conj.state
                match = re.match(r'\s*(\d+)\s+(\d+)\s+([-\d.]+)\s+([-\d]+)\s+([-\d]+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([\dE.+-]+)\s+([\d.]+)\s+(\d+)', line)
                
                if match:
                    number = int(match.group(1))
                    isospin = int(match.group(2))  # 1=proton, 2=neutron
                    jz = float(match.group(3))
                    K = int(match.group(4))  # tot quantum number
                    v2 = float(match.group(8))
                    energy = float(match.group(11))
                    conjugate = int(match.group(14))
                    
                    states.append({
                        'number': number,
                        'isospin': isospin,
                        'jz': jz,
                        'K': K,
                        'v2': v2,
                        'energy': energy,
                        'conjugate': conjugate
                    })
        
        return states
    
    def find_candidate_states(self, states, isospin, max_states=10):
        """
        Find candidate blocking states from spectrum.
        
        Filters conjugate pairs (prefers positive j_z when possible).
        
        Args:
            states: List of state dictionaries
            isospin: 1 for proton, 2 for neutron
            max_states: Maximum number of candidates to return
            
        Returns:
            list: Candidate states sorted by proximity to Fermi surface
        """
        candidates = []
        seen_pairs = set()  # Track conjugate pairs to avoid duplicates
        
        for state in states:
            if state['isospin'] != isospin:
                continue
            
            # Skip if we've already added the conjugate partner
            pair = tuple(sorted([state['number'], state['conjugate']]))
            if pair in seen_pairs:
                continue
            
            # For blocking, look for states with significant occupation
            # v² close to 0.5 indicates Fermi surface
            if 0.1 < state['v2'] < 0.99:
                # Prefer positive j_z from conjugate pair
                # Check if conjugate exists and has opposite j_z
                conj_state = next((s for s in states if s['number'] == state['conjugate']), None)
                
                if conj_state and conj_state['jz'] * state['jz'] < 0:
                    # Conjugate pair with opposite j_z - prefer positive
                    if state['jz'] < 0:
                        continue  # Skip this one, the positive j_z version will be added
                
                # Score by distance from v²=0.5 (Fermi surface)
                score = abs(state['v2'] - 0.5)
                candidates.append((score, state))
                seen_pairs.add(pair)
        
        # Sort by score (closest to Fermi surface first)
        candidates.sort(key=lambda x: x[0])
        
        # Return just the state dictionaries, not scores
        return [state for score, state in candidates[:max_states]]
    
    @staticmethod
    def test_single_blocking(args):
        """
        Static method to test a single blocking configuration.
        Used for parallel execution.
        
        Args:
            args: Tuple of (scanner, state, isospin, i, total, ref_energy)
        
        Returns:
            dict: Results including energy and state info
        """
        scanner, state, isospin, i, total, ref_energy = args
        
        istblo_p = state['number'] if isospin == 1 else 0
        istblo_n = state['number'] if isospin == 2 else 0
        
        iso_label = 'p' if isospin == 1 else 'n'
        log_prefix = f"block_{state['number']}"
        
        print(f"  [{i}/{total}] Testing istblo=({istblo_p},{istblo_n}) "
              f"[{iso_label}, K={state['K']}, v²={state['v2']:.3f}, jz={state['jz']:.1f}]")
        
        # Create input with blocking (use unique subdirectory to avoid conflicts)
        state_dir = scanner.workdir / f"state_{state['number']}"
        state_dir.mkdir(exist_ok=True, parents=True)
        test_input = state_dir / "for005"
        scanner.modify_input(scanner.input_file, test_input,
                           istblo_p=istblo_p, istblo_n=istblo_n)
        
        # Temporarily change workdir for this calculation
        original_workdir = scanner.workdir
        scanner.workdir = state_dir
        energy = scanner.run_skyax(test_input, log_prefix)
        scanner.workdir = original_workdir
        
        # Move results back to main workdir
        if energy is not None:
            for f in state_dir.glob(f"{log_prefix}_*"):
                dest = scanner.workdir / f.name
                if dest.exists():
                    dest.unlink()
                shutil.move(str(f), str(dest))
        
        result = {
            'state_number': state['number'],
            'isospin': isospin,
            'K': state['K'],
            'jz': state['jz'],
            'v2': state['v2'],
            'istblo_p': istblo_p,
            'istblo_n': istblo_n,
            'energy': energy
        }
        
        if energy is not None:
            delta = energy - ref_energy
            print(f"    → Energy: {energy:.3f} MeV (ΔE = {delta:+.3f} vs no blocking)")
        else:
            print(f"    → FAILED")
        
        return result
    
    def scan_blocking_states(self, max_states=10):
        """
        Main workflow with parallel state testing.
        
        Args:
            max_states: Maximum number of candidate states to test
            
        Returns:
            dict: Best configuration with istblo_p, istblo_n, energy
        """
        print("=" * 70)
        print("AUTOMATED BLOCKING STATE SELECTION")
        print("=" * 70)
        print(f"Parallelization: {self.num_processes} Python processes × "
              f"{self.num_threads} OpenMP threads = {self.num_processes * self.num_threads} total threads")
        
        # Parse target nucleus
        target_nprot, target_nneut = self.parse_input_nucleus(self.input_file)
        
        odd_protons = (target_nprot % 2 == 1)
        odd_neutrons = (target_nneut % 2 == 1)
        
        print(f"\nTarget nucleus: Z={target_nprot}, N={target_nneut}")
        print(f"  Configuration: {'odd-Z' if odd_protons else 'even-Z'}, "
              f"{'odd-N' if odd_neutrons else 'even-N'}")
        
        if not (odd_protons or odd_neutrons):
            print("\nWARNING: Even-even nucleus detected. No blocking needed.")
            return {'istblo_p': 0, 'istblo_n': 0, 'energy': None}
        
        # Step 1: Run target nucleus with NO blocking
        print("\n" + "=" * 70)
        print("STEP 1: Initial Calculation (No Blocking)")
        print("=" * 70)
        print(f"Running {target_nprot}p + {target_nneut}n with istblo=(0,0)...")
        
        ref_input = self.workdir / "for005_initial"
        self.modify_input(self.input_file, ref_input, istblo_p=0, istblo_n=0)
        
        print("  Running: initial...")
        ref_energy = self.run_skyax(ref_input, "initial")
        
        if ref_energy is None:
            print("\nERROR: Initial calculation failed")
            return None
        
        print(f"  Initial energy (no blocking): {ref_energy:.3f} MeV")
        
        # Parse spectrum
        details_file = self.workdir / "initial_out_details.sk2"
        if not details_file.exists():
            print(f"\nERROR: Details file not found: {details_file}")
            return None
        
        states = self.parse_state_details(details_file)
        print(f"  Parsed {len(states)} states from spectrum")
        
        # Step 2: Find candidates
        print("\n" + "=" * 70)
        print("STEP 2: Identify Blocking Candidates")
        print("=" * 70)
        
        candidates = []
        
        if odd_protons:
            p_cand = self.find_candidate_states(states, isospin=1, max_states=max_states)
            candidates.extend([(state, 1) for state in p_cand])
            print(f"\nProton candidates (states near Fermi surface, positive j_z preferred):")
            for i, state in enumerate(p_cand, 1):
                print(f"  {i}. State {state['number']:3d}: K={state['K']:2d}, jz={state['jz']:4.1f}, "
                      f"v²={state['v2']:.3f}, E={state['energy']:7.2f} MeV, conj={state['conjugate']:3d}")
        
        if odd_neutrons:
            n_cand = self.find_candidate_states(states, isospin=2, max_states=max_states)
            candidates.extend([(state, 2) for state in n_cand])
            print(f"\nNeutron candidates (states near Fermi surface, positive j_z preferred):")
            for i, state in enumerate(n_cand, 1):
                print(f"  {i}. State {state['number']:3d}: K={state['K']:2d}, jz={state['jz']:4.1f}, "
                      f"v²={state['v2']:.3f}, E={state['energy']:7.2f} MeV, conj={state['conjugate']:3d}")
        
        if not candidates:
            print("\nWARNING: No blocking candidates found!")
            return {'istblo_p': 0, 'istblo_n': 0, 'energy': None}
        
        print(f"\nTotal configurations to test: {len(candidates)}")
        
        # Step 3: Test blocking configurations (in parallel)
        print("\n" + "=" * 70)
        print(f"STEP 3: Test Blocking Configurations (Parallel)")
        print("=" * 70)
        
        # Prepare arguments for parallel execution
        test_args = [
            (self, state, isospin, i, len(candidates), ref_energy)
            for i, (state, isospin) in enumerate(candidates, 1)
        ]
        
        # Run in parallel
        if self.num_processes > 1:
            with Pool(processes=self.num_processes) as pool:
                results = pool.map(BlockingScanner.test_single_blocking, test_args)
        else:
            # Serial execution
            results = [BlockingScanner.test_single_blocking(args) for args in test_args]
        
        # Filter out failed calculations
        results = [r for r in results if r['energy'] is not None]
        
        if not results:
            print("\nERROR: All blocking calculations failed")
            return None
        
        # Step 4: Results
        print("\n" + "=" * 70)
        print("RESULTS SUMMARY")
        print("=" * 70)
        
        print(f"\nInitial (no blocking): E = {ref_energy:.3f} MeV")
        
        best = min(results, key=lambda x: x['energy'])
        
        print("\nBlocking configurations tested (sorted by energy):")
        for r in sorted(results, key=lambda x: x['energy']):
            delta = r['energy'] - best['energy']
            delta_vs_initial = r['energy'] - ref_energy
            marker = " <-- BEST" if r == best else ""
            iso_label = 'p' if r['isospin'] == 1 else 'n'
            print(f"  istblo=({r['istblo_p']:3d},{r['istblo_n']:3d}) [{iso_label}, K={r['K']:2d}, jz={r['jz']:4.1f}]: "
                  f"E = {r['energy']:10.3f} MeV  (ΔE = {delta:+7.3f} MeV, "
                  f"vs no-block: {delta_vs_initial:+7.3f}){marker}")
        
        print(f"\n{'=' * 70}")
        print(f"RECOMMENDED: istblo=({best['istblo_p']},{best['istblo_n']})")
        print(f"  State {best['state_number']}: K={best['K']}, jz={best['jz']:.1f}, v²={best['v2']:.3f}")
        print(f"Ground state energy: {best['energy']:.3f} MeV")
        print(f"Improvement over no blocking: {best['energy'] - ref_energy:+.3f} MeV")
        print(f"{'=' * 70}")
        
        return best


def main():
    parser = argparse.ArgumentParser(
        description='Automated blocking state selection for SkyAx odd nuclei',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Physical approach:
  - Runs target odd nucleus with NO blocking to identify spectrum
  - Finds states near Fermi surface (partially occupied, v² ~ 0.5)
  - Tests blocking each candidate state (in parallel)
  - Selects configuration with lowest energy

Parallelization:
  --processes: Number of parallel Python processes (tests multiple states simultaneously)
  --threads: Number of OpenMP threads per process
  Total parallelism = processes × threads

Example:
  # 2 Python processes, 2 OpenMP threads each = 4 total threads
  python auto_blocking_v3.py --input for005 --exec ../../skyax.gfortran.omp \\
                              --processes 2 --threads 2 --max-states 6
        """
    )
    
    parser.add_argument('--input', '-i', required=True,
                        help='Input file for target odd nucleus')
    parser.add_argument('--exec', '-e', required=True,
                        help='Path to SkyAx executable')
    parser.add_argument('--workdir', '-w', default='autoblock_v3',
                        help='Working directory (default: autoblock_v3)')
    parser.add_argument('--threads', '-t', type=int, default=1,
                        help='Number of OpenMP threads per process (default: 1)')
    parser.add_argument('--processes', '-p', type=int, default=1,
                        help='Number of parallel Python processes (default: 1)')
    parser.add_argument('--max-states', '-m', type=int, default=10,
                        help='Maximum candidate states to test (default: 10)')
    
    args = parser.parse_args()
    
    # Validate parallelization
    total_threads = args.processes * args.threads
    available_cpus = cpu_count()
    if total_threads > available_cpus:
        print(f"WARNING: Total threads ({total_threads}) exceeds available CPUs ({available_cpus})")
        print(f"         Consider reducing --processes or --threads")
    
    scanner = BlockingScanner(
        args.input,
        args.exec,
        args.workdir,
        args.threads,
        args.processes
    )
    
    result = scanner.scan_blocking_states(max_states=args.max_states)
    
    if result and result['energy'] is not None:
        # Save best configuration
        final_input = Path(args.workdir) / "for005_best"
        scanner.modify_input(
            args.input,
            final_input,
            istblo_p=result['istblo_p'],
            istblo_n=result['istblo_n']
        )
        print(f"\nBest configuration saved to: {final_input}")
        return 0
    else:
        print("\nERROR: Blocking scan failed")
        return 1


if __name__ == '__main__':
    sys.exit(main())
