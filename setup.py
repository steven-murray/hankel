from setuptools import setup
import os
import sys
import re
import io

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

if sys.argv[-1] == "publish":
    os.system("rm dist/*")
    os.system("python setup.py sdist")
    os.system("python setup.py bdist_wheel")
    os.system("twine upload dist/*")
    sys.exit()

setup(name='hankel',
      install_requires=['numpy>=1.6.1', 'scipy>=0.12.0', 'mpmath>=0.18'],
      version = find_version("hankel", "__init__.py"),
      packages=['hankel'],
      description='Hankel Transformations using method of Ogata 2005',
      long_description=read('README.rst'),
      author='Steven Murray',
      author_email='steven.murray@curtin.edu.au',
      license = "MIT",
      url='https://github.com/steven-murray/hankel')
