from setuptools import setup
import os
import sys

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    os.system("python setup.py bdist_wheel upload")
    sys.exit()

setup(name='hankel',
      install_requires=['numpy', 'scipy', 'mpmath'],
      version='0.2.2',
      py_modules=['hankel'],
      description='Hankel Transformations using method of Ogata 2005',
      long_description=read('README.rst'),
      author='Steven Murray',
      author_email='steven.murray@curtin.edu.au',
      url='https://github.com/steven-murray/hankel')
