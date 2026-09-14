"""Run HMMER executables."""

from __future__ import annotations

import logging
import shutil
import subprocess

from biocol.exceptions import HmmExecutionError

logger = logging.getLogger(__name__)


def run_hmmer(command: list[str]) -> None:
    executable = command[0]
    if shutil.which(executable) is None:
        raise HmmExecutionError(
            f"'{executable}' was not found in PATH. Install HMMER in the conda environment."
        )
    logger.debug("HMMER command: %s", " ".join(command))
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise HmmExecutionError(
            f"{executable} failed (exit code {completed.returncode}): {completed.stderr.strip()}"
        )
