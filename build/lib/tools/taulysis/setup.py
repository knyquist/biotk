import sys
import os.path
from setuptools import setup, Extension, find_packages

REQUIREMENTS_TXT = "requirements.txt"

if("install" in sys.argv) and sys.version_info < (2, 7, 0):
    print "taulysis requires Python 2.7"
    sys.exit(-1)

globals = {}
execfile("taulysis/__version__.py", globals)
__VERSION__ = globals["__VERSION__"]

setup(name="taulysis",
      version=__VERSION__,
      author="Kristofor Nyquist",
      author_email="knyquist@pacificbiosciences.com",
      description="Estimate first-pass, second-pass, and rolling-circle tau",
      packages=find_packages('.'),
      zip_safe=True,
      entry_points={
          "console_scripts" : ["taulysis = taulysis.main:main"]})#,
      # install_requires=_get_local_requirements(REQUIREMENTS_TXT))
