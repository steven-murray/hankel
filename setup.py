from setuptools import setup
import os
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(name='hankel',
      install_requires=['numpy', 'scipy', 'mpmath'],
      version='0.2.0',
      py_modules=['hankel'],
      description='Hankel Transformations using method of Ogata 2005',
      long_description=read('readme.rst'),
      author='Steven Murray',
      author_email='steven.murray@uwa.edu.au',
      url='https://github.com/steven-murray/hankel')
