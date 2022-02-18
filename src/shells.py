#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from subprocess import CalledProcessError
import os, subprocess

# Local imports
from utils import fatal, err


def set_options(strict):
    """
    Changes behavior of default shell and get overrides options 
    to run bash in a strict mode. 
    @param strict <bool>:
        Overrides default shell options and runs shell in strict or 
        less permissive mode.
    @return prefix <int>:
        Returns overrides options to run bash in a strict mode
    """
    prefix = ''  # permissive shell option
    if strict: 
        # Changes behavior of default shell
        # set -e: exit immediately upon error
        # set -u: treats unset variables as an error
        # set -o pipefail: exits if a error occurs in any point of a pipeline
        prefix = 'set -euo pipefail; '

    return prefix


def bash(cmd, interpreter='/bin/bash', strict=set_options(True), cwd=os.getcwd(), **kwargs):
    """
    Interface to run a process or bash command. Using subprocess.call_check()
    due to portability across most python versions. It was introduced in python 2.5
    and it is also interoperabie across all python 3 versions. 
    @param cmd <str>:
        Shell command to run
    @param interpreter <str>:
        Interpreter for command to run [default: bash]
    @pararm strict <bool>:
        Prefixes any command with 'set -euo pipefail' to ensure process fail with
        the expected exit-code  
    @params kwargs <check_call()>:
        Keyword arguments to modify subprocess.check_call() behavior
    @return exitcode <int>:
        Returns the exit code of the run command, failures return non-zero exit codes
    """
    try:
        exitcode = subprocess.check_call(strict + cmd, 
            shell=True, 
            executable=interpreter, 
            cwd=cwd, 
            **kwargs
        )
    except CalledProcessError as e:
        exitcode = e.returncode
        err("""WARNING: Failed to run '{}' command!
        └── Command returned a non-zero exitcode of '{}'.""".format(strict + cmd, exitcode)
        )

    return exitcode


if __name__ == '__main__':
    # Tests
    bash('ls -la /home/')
    bash('ls -la /fake/dne/path')
