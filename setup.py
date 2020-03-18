"""hankel: Hankel Transformations using method of Ogata 2005."""

import io
import os
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ) as fp:
        return fp.read()


CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Developers",
    "Intended Audience :: End Users/Desktop",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Natural Language :: English",
    "Operating System :: Unix",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3 :: Only",
    "Topic :: Scientific/Engineering",
    "Topic :: Utilities",
]

with open("requirements_dev.txt") as fl:
    test_req = fl.readlines()

with open("docs/requirements.txt") as fl:
    doc_req = fl.readlines()

# use first version of req. to support py35
setup(
    name="hankel",
    install_requires=["numpy>=1.9.3", "scipy>=0.16.1", "mpmath>=1.0.0"],
    classifiers=CLASSIFIERS,
    python_requires=">=3.5",
    packages=["hankel"],
    description="Hankel Transformations using method of Ogata 2005",
    long_description=read("README.rst"),
    author="Steven Murray",
    author_email="steven.murray@curtin.edu.au",
    license="MIT",
    extras_require={"dev": test_req + doc_req, "tests": test_req, "docs": doc_req},
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    url="https://github.com/steven-murray/hankel",
)
