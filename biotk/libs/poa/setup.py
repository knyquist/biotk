from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'Levenshtein Distance',
  ext_modules=[
    Extension('levenshtein_distance',
              sources=['levenshtein_distance.pyx'],
              extra_compile_args=['-fgnu89-inline'],
              language='c')
    ],
  cmdclass = {'build_ext': build_ext}
)