#!/usr/bin/env python3
"""
Run all COBRAxy tests.

This script runs all test files and provides a summary of results.
Can be run with or without pytest.
"""

import sys
import os
import subprocess
from pathlib import Path

# Colors for terminal output
class Colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def print_header(text):
    """Print a formatted header"""
    print(f"\n{Colors.HEADER}{Colors.BOLD}{'='*70}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{text:^70}{Colors.ENDC}")
    print(f"{Colors.HEADER}{Colors.BOLD}{'='*70}{Colors.ENDC}\n")


def run_with_pytest():
    """Run tests using pytest"""
    print_header("Running Tests with pytest")
    
    test_dir = Path(__file__).parent
    
    # Run pytest with coverage
    cmd = [
        sys.executable, "-m", "pytest",
        str(test_dir),
        "-v",
        "--tb=short",
        f"--cov={test_dir.parent}",
        "--cov-report=term-missing",
        "--cov-report=html"
    ]
    
    print(f"{Colors.OKCYAN}Command: {' '.join(cmd)}{Colors.ENDC}\n")
    
    result = subprocess.run(cmd)
    
    if result.returncode == 0:
        print(f"\n{Colors.OKGREEN}{Colors.BOLD}✓ All tests passed!{Colors.ENDC}")
        print(f"{Colors.OKGREEN}Coverage report generated in htmlcov/index.html{Colors.ENDC}")
    else:
        print(f"\n{Colors.FAIL}{Colors.BOLD}✗ Some tests failed!{Colors.ENDC}")
    
    return result.returncode


def run_without_pytest():
    """Run tests without pytest"""
    print_header("Running Tests (without pytest)")
    print(f"{Colors.WARNING}Note: pytest not found. Running basic test execution.{Colors.ENDC}\n")
    
    test_dir = Path(__file__).parent
    test_files = [
        "test_utils.py",
        "test_generators.py",
        "test_marea.py",
        "test_clustering.py"
    ]
    
    total_passed = 0
    total_failed = 0
    failed_files = []
    
    for test_file in test_files:
        test_path = test_dir / test_file
        
        if not test_path.exists():
            print(f"{Colors.WARNING}⊘ {test_file}: Not found{Colors.ENDC}")
            continue
        
        print(f"\n{Colors.OKBLUE}{Colors.BOLD}Running {test_file}...{Colors.ENDC}")
        print(f"{Colors.OKBLUE}{'─'*70}{Colors.ENDC}")
        
        result = subprocess.run(
            [sys.executable, str(test_path)],
            capture_output=False
        )
        
        if result.returncode == 0:
            print(f"{Colors.OKGREEN}✓ {test_file} passed{Colors.ENDC}")
        else:
            print(f"{Colors.FAIL}✗ {test_file} failed{Colors.ENDC}")
            failed_files.append(test_file)
    
    # Print summary
    print_header("Test Summary")
    
    if not failed_files:
        print(f"{Colors.OKGREEN}{Colors.BOLD}✓ All test files passed!{Colors.ENDC}")
        return 0
    else:
        print(f"{Colors.FAIL}{Colors.BOLD}✗ Failed test files:{Colors.ENDC}")
        for file in failed_files:
            print(f"{Colors.FAIL}  - {file}{Colors.ENDC}")
        return 1


def check_dependencies():
    """Check if required dependencies are installed"""
    print_header("Checking Dependencies")
    
    required = [
        "pandas", "numpy", "scipy", "sklearn", 
        "cobra", "lxml", "matplotlib", "seaborn"
    ]
    
    missing = []
    
    for package in required:
        try:
            __import__(package)
            print(f"{Colors.OKGREEN}✓ {package:20} installed{Colors.ENDC}")
        except ImportError:
            print(f"{Colors.FAIL}✗ {package:20} missing{Colors.ENDC}")
            missing.append(package)
    
    if missing:
        print(f"\n{Colors.WARNING}Missing packages: {', '.join(missing)}{Colors.ENDC}")
        print(f"{Colors.WARNING}Install with: pip install {' '.join(missing)}{Colors.ENDC}")
        return False
    
    return True


def main():
    """Main entry point"""
    print_header("COBRAxy Test Suite")
    
    # Check dependencies
    if not check_dependencies():
        print(f"\n{Colors.FAIL}Please install missing dependencies first.{Colors.ENDC}")
        return 1
    
    # Try to use pytest if available
    try:
        import pytest
        return run_with_pytest()
    except ImportError:
        return run_without_pytest()


if __name__ == "__main__":
    sys.exit(main())
