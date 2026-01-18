#!/usr/bin/env python3
"""
Automated blocking state selection for SkyAx odd-nucleus calculations.

Physical principle:
- Run the target odd nucleus with NO blocking (istblo=0,0)
- Identify states near the Fermi surface (partial occupation v²~0.5)
- Test blocking each candidate state
- Select configuration with lowest energy

The blocked state prevents pairing correlations for that orbital,
representing the unpaired nucleon in the odd system.
"""

import sys
import re
import subprocess
import shutil
from pathlib import Path
import argparse


class BlockingScanner:
    def __init__(self, input_file, executable, workdir='.', num_threads=1):
        """
        Initialize the blocking scanner.
        
        Args:
            input_file: Path to input file (for TARGET odd nucleus)
            executable: Path to SkyAx executable
            workdir: Working directory for calculations
            num_threads: Number of OpenMP threads
        """
        self.input_file = Path(input_file).resolve()
        self.executable = Path(executable).resolve()
        self.workdir = Path(workdir).resolve()
        self.num_threads = num_threads
        
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
    
    def modify_input(self, input_file, output_file, istblo_p=None, istblo_n=None, 
                    nprot=None, nneut=None):
        """
        Modify input file to set blocking states and/or particle numbers.
        
        Args:
            input_file: Original input file path
            output_file: Path for modified input file
            istblo_p: Proton blocking state number (None to keep original)
            istblo_n: Neutron blocking state number (None to keep original)
            nprot: Proton number (None to keep original)
            nneut: Neutron number (None to keep original)
        """
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Replace particle numbers if specified
        if nprot is not None:
            content = re.sub(r'nprot\s*=\s*\d+', f'nprot={nprot}', content)
        if nneut is not None:
            content = re.sub(r'nneut\s*=\s*\d+', f'nneut={nneut}', content)
        
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
        import os
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(self.num_threads)
        
        log_file = self.workdir / f"{log_prefix}.log"
        
        print(f"  Running: {log_file.stem}...")
        
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
            
            # Archive output files
            for pattern in ['out_*']:
                for f in self.workdir.glob(pattern):
                    new_name = f.with_name(f"{log_prefix}_{f.name}")
                    if new_name.exists():
                        new_name.unlink()
                    shutil.move(str(f), str(new_name))
            
            return energy
            
        except subprocess.TimeoutExpired:
            print(f"    WARNING: Run timed out")
            return None
        except Exception as e:
            print(f"    ERROR: {e}")
            return None
    
    def parse_state_details(self, details_file):
        """
        Parse out_details.sk2 file to extract single-particle state information.
        
        Returns:
            list: List of state dictionaries with keys:
                  number, isospin (iq), jz, K, v2, energy, conjugate, etc.
        """
        states = []
        
        with open(details_file, 'r') as f:
            for line in f:
                # Format: "  Nr Iso j_z tot r z wtpair v² uv delta energy ... conj"
                # Example: "  398 2  1.5   4   1   4   0.99997 0.83282 ... 290"
                # We need to capture more fields to get the conjugate (last column)
                match = re.match(r'\s*(\d+)\s+(\d+)\s+([-\d.]+)\s+([-\d]+)\s+\d+\s+([-\d]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*\s+(\d+)\s*$', line)
                if match:
                    number = int(match.group(1))
                    isospin = int(match.group(2))  # 1=proton, 2=neutron
                    jz = float(match.group(3))
                    K = int(match.group(4))
                    v2 = float(match.group(7))
                    energy = float(match.group(10))
                    conjugate = int(match.group(11))
                    
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
        
        Look for states near Fermi surface with partial occupation (v² ~ 0.5)
        or highly occupied states (v² > 0.5) that represent the unpaired nucleon.
        
        Filters out conjugate pairs (only keeps one state from each conjugate pair).
        
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
            # v² > 0.5 indicates occupied states that could be blocked
            if 0.1 < state['v2'] < 0.99:
                # Score by distance from v²=0.5 (Fermi surface)
                score = abs(state['v2'] - 0.5)
                candidates.append((score, state))
                seen_pairs.add(pair)
        
        # Sort by score (closest to Fermi surface first)
        candidates.sort(key=lambda x: x[0])
        
        # Return just the state dictionaries, not scores
        return [state for score, state in candidates[:max_states]]
    
    def scan_blocking_states(self, max_states=10):
        """
        Main workflow:
        1. Determine target nucleus (odd-N, odd-Z, or odd-odd)
        2. Run target nucleus with NO blocking to get reference spectrum
        3. Extract candidate states (near Fermi surface) from this spectrum
        4. Test each blocking configuration
        5. Return best configuration
        
        Args:
            max_states: Maximum number of candidate states to test
            
        Returns:
            dict: Best configuration with istblo_p, istblo_n, energy
        """
        print("=" * 70)
        print("AUTOMATED BLOCKING STATE SELECTION")
        print("=" * 70)
        
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
        
        # Step 1: Run target nucleus with NO blocking to get reference spectrum
        print("\n" + "=" * 70)
        print("STEP 1: Initial Calculation (No Blocking)")
        print("=" * 70)
        print(f"Running {target_nprot}p + {target_nneut}n with istblo=(0,0) to identify candidates...")
        
        ref_input = self.workdir / "for005_initial"
        self.modify_input(self.input_file, ref_input, 
                         istblo_p=0, istblo_n=0)
        
        ref_energy = self.run_skyax(ref_input, "initial")
        
        if ref_energy is None:
            print("\nERROR: Initial calculation failed")
            return None
        
        print(f"  Initial energy (no blocking): {ref_energy:.3f} MeV")
        
        # Parse spectrum from unblocked calculation
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
            print(f"\nProton candidates (states near Fermi surface):")
            for i, state in enumerate(p_cand, 1):
                print(f"  {i}. State {state['number']:3d}: K={state['K']:2d}, "
                      f"v²={state['v2']:.3f}, E={state['energy']:7.2f} MeV, "
                      f"conj={state['conjugate']:3d}")
        
        if odd_neutrons:
            n_cand = self.find_candidate_states(states, isospin=2, max_states=max_states)
            candidates.extend([(state, 2) for state in n_cand])
            print(f"\nNeutron candidates (states near Fermi surface):")
            for i, state in enumerate(n_cand, 1):
                print(f"  {i}. State {state['number']:3d}: K={state['K']:2d}, "
                      f"v²={state['v2']:.3f}, E={state['energy']:7.2f} MeV, "
                      f"conj={state['conjugate']:3d}")
        
        if not candidates:
            print("\nWARNING: No blocking candidates found!")
            return {'istblo_p': 0, 'istblo_n': 0, 'energy': None}
        
        print(f"\nTotal configurations to test: {len(candidates)}")
        
        # Step 3: Test blocking configurations
        print("\n" + "=" * 70)
        print(f"STEP 3: Test Blocking Configurations")
        print("=" * 70)
        
        results = []
        
        for i, (state, isospin) in enumerate(candidates, 1):
            istblo_p = state['number'] if isospin == 1 else 0
            istblo_n = state['number'] if isospin == 2 else 0
            
            iso_label = 'p' if isospin == 1 else 'n'
            print(f"\nConfig {i}/{len(candidates)}: istblo=({istblo_p},{istblo_n}) "
                  f"[{iso_label}, K={state['K']}, v²={state['v2']:.3f}]")
            
            # Create input with blocking
            test_input = self.workdir / f"for005_block_{state['number']}"
            self.modify_input(self.input_file, test_input,
                            istblo_p=istblo_p, istblo_n=istblo_n)
            
            energy = self.run_skyax(test_input, f"block_{state['number']}")
            
            if energy is not None:
                results.append({
                    'state_number': state['number'],
                    'isospin': isospin,
                    'K': state['K'],
                    'v2': state['v2'],
                    'istblo_p': istblo_p,
                    'istblo_n': istblo_n,
                    'energy': energy
                })
                delta = energy - ref_energy
                print(f"  Energy: {energy:.3f} MeV (ΔE = {delta:+.3f} vs no blocking)")
            else:
                print(f"  FAILED")
        
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
            print(f"  istblo=({r['istblo_p']:3d},{r['istblo_n']:3d}) [{iso_label}]: "
                  f"E = {r['energy']:10.3f} MeV  (ΔE = {delta:+7.3f} MeV, "
                  f"vs no-block: {delta_vs_initial:+7.3f}){marker}")
        
        print(f"\n{'=' * 70}")
        print(f"RECOMMENDED: istblo=({best['istblo_p']},{best['istblo_n']})")
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
  - Tests blocking each candidate state
  - Selects configuration with lowest energy

Example:
  python auto_blocking_v2.py --input for005 --exec ../../skyax.gfortran.omp --threads 4 --max-states 5
        """
    )
    
    parser.add_argument('--input', '-i', required=True,
                        help='Input file for target odd nucleus')
    parser.add_argument('--exec', '-e', required=True,
                        help='Path to SkyAx executable')
    parser.add_argument('--workdir', '-w', default='autoblock_v2',
                        help='Working directory (default: autoblock_v2)')
    parser.add_argument('--threads', '-t', type=int, default=1,
                        help='Number of OpenMP threads (default: 1)')
    parser.add_argument('--max-states', '-m', type=int, default=10,
                        help='Maximum candidate states to test (default: 10)')
    
    args = parser.parse_args()
    
    scanner = BlockingScanner(
        args.input,
        args.exec,
        args.workdir,
        args.threads
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
