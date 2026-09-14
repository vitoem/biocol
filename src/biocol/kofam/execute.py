"""Run the KofamScan ``exec_annotation`` executable."""

from __future__ import annotations

import logging
import shutil
import subprocess

from biocol.exceptions import KofamExecutionError

logger = logging.getLogger(__name__)


def run_exec_annotation(command: list[str]) -> None:
    executable = command[0]
    if shutil.which(executable) is None:
        raise KofamExecutionError(
            f"'{executable}' was not found in PATH. Install KofamScan "
            "(exec_annotation) in the conda environment."
        )
    logger.debug("KofamScan command: %s", " ".join(command))
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise KofamExecutionError(
            f"{executable} failed (exit code {completed.returncode}): {completed.stderr.strip()}"
        )
