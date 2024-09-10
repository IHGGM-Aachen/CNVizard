"""
File to execute the CNVizard
@author: Jeremias Krause , Carlos Classen, Matthias Begemann. Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import sys
import argparse
import subprocess
from pathlib import Path

def main():
    """
    Entry point for launching the CNVizard Streamlit app.
    """
    # Parse command-line arguments, allowing an optional .env file
    parser = argparse.ArgumentParser(description="Run CNVizard Streamlit app.")
    parser.add_argument("env", nargs="?", default=None, help="Path to the .env file.")  # Optional .env file
    args = parser.parse_args()

    # Get the absolute path to the app.py (this script is assumed to be located in the same directory as app.py)
    app_path = Path(__file__).parent / "app.py"

    # Get the Python executable that is currently running this script
    executable = sys.executable

    # Build the command to launch the Streamlit app using the correct Python executable
    command = [str(executable), "-m", "streamlit", "run", str(app_path)]

    # If an environment file is provided, add it as an argument to the command
    if args.env:
        command.append(str(args.env))

    # Run the command to launch Streamlit as a subprocess
    # The 'check=True' ensures that if the subprocess fails, it raises an exception
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
