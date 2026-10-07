import pathlib
import subprocess
import sys

from setuptools import setup
from setuptools.command.build import build as _build
from setuptools.command.build_ext import build_ext as _build_ext


ROOT = pathlib.Path(__file__).resolve().parent


def run_make():
    make_cmd = 'make'
    if sys.platform.startswith('win'):
        make_cmd = 'mingw32-make'
    subprocess.check_call([make_cmd], cwd=str(ROOT))


class BuildWithMake(_build):
    def run(self):
        run_make()
        super().run()


class BuildExtWithMake(_build_ext):
    def run(self):
        run_make()


setup(
    name='spinOS-kepler',
    version='0.1.0',
    description='Build helper for spinOS C Kepler solver',
    packages=[],
    py_modules=[],
    cmdclass={
        'build': BuildWithMake,
        'build_ext': BuildExtWithMake,
    },
)
