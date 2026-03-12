#!/usr/bin/env python3
"""
Simple wrapper for auto_blocking_v3.py with sensible defaults.
"""

import sys
import argparse
from pathlib import Path

# Import the main module
sys.path.insert(0, str(Path(__file__).parent))
from auto_blocking_v3 import BlockingScanner


def main():
    parser = argparse.ArgumentParser(
        description='Quick blocking scan with sensible defaults'
    )
    
    parser.add_argument('input_file', help='Input file')
    parser.add_argument('--exec', '-e', default='../../skyax.gfortran.omp',
                        help='Executable (default: ../../skyax.gfortran.omp)')
    parser.add_argument('--workdir', '-w', default='autoblock',
                        help='Working directory (default: autoblock)')
    parser.add_argument('--cpus', '-c', type=int, default=4,
                        help='Total CPUs to use (default: 4)')
    parser.add_argument('--max-states', '-m', type=int, default=6,
                        help='Max states to test (default: 6)')
    parser.add_argument('--test', action='store_true',
                        help='Quick test mode (2 states only)')
    
    args = parser.parse_args()
    
    # Set parallelization based on total CPUs
    if args.test:
        processes = 2
        threads = 1
        max_states = 2
    else:
        # Use balanced parallelization
        if args.cpus >= 8:
            processes = 4
            threads = args.cpus // 4
        elif args.cpus >= 4:
            processes = 2
            threads = args.cpus // 2
        else:
            processes = 1
            threads = args.cpus
        max_states = args.max_states
    
    print(f"Configuration: {processes} processes × {threads} threads = {processes*threads} CPUs")
    
    scanner = BlockingScanner(
        args.input_file,
        args.exec,
        args.workdir,
        threads,
        processes
    )
    
    result = scanner.scan_blocking_states(max_states=max_states)
    
    if result and result['energy'] is not None:
        final_input = Path(args.workdir) / "for005_best"
        scanner.modify_input(
            args.input_file,
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
