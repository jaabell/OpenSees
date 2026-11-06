"""Compile and run tests/test_logstrain_kernel.cpp (also registered with CTest).

The C++ driver checks SRC/material/nD/LogStrainKernel.h, the numerical core of
nDMaterial LogStrain, against closed-form Hencky stresses, frame indifference
and a central-difference spatial modulus at distinct and repeated principal
stretches. The header has no OpenSees dependency, so the test needs only a C++
compiler (g++, clang++ or c++ on PATH, or MSVC cl in a developer shell) and is
skipped when none is available.
"""
import os
import shutil
import subprocess

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "test_logstrain_kernel.cpp")
_INC = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))


def _compile(tmp_path):
    exe = str(tmp_path / ("logstrain_kernel" + (".exe" if os.name == "nt" else "")))
    for cxx in ("g++", "clang++", "c++"):
        path = shutil.which(cxx)
        if path:
            cmd = [path, "-O2", "-std=c++11", "-I", _INC, _SRC, "-o", exe]
            break
    else:
        path = shutil.which("cl")
        if not path:
            pytest.skip("no C++ compiler on PATH")
        cmd = [path, "/nologo", "/EHsc", "/O2", "/I", _INC, _SRC,
               "/Fe:" + exe, "/Fo:" + str(tmp_path) + os.sep]
    cc = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True,
                        text=True, cwd=str(tmp_path))
    assert cc.returncode == 0, f"compilation failed:\n{cc.stdout}\n{cc.stderr}"
    return exe


def test_logstrain_kernel(tmp_path):
    exe = _compile(tmp_path)
    run = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True,
                         text=True, timeout=60)
    assert run.returncode == 0, run.stdout + run.stderr
    assert "all LogStrainKernel checks passed" in run.stdout
