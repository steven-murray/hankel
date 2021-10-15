# .. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs by opening a `new issue <https://github.com/steven-murray/hankel/issues/new?template=bugs.md>`_.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with
`"bug" and "help wanted" <https://github.com/steven-murray/hankel/issues?q=is%3Aissue+is%3Aopen+label%3Abug+label%3A%22help+wanted%22+>`_
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with
`"enhancement" and "help wanted" <https://github.com/steven-murray/hankel/issues?utf8=%E2%9C%93&q=is%3Aissue+is%3Aopen+label%3Aenhancement+label%3A%22help+wanted%22+>`_
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

hankel could always use more documentation, whether as part of the
official hankel docs, in docstrings, or even on the web in blog posts,
articles, and such.

In particular, if you have used hankel in a novel way, please consider submitting a
script or Jupyter notebook as "demo" documentation.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to
`file an issue <https://github.com/steven-murray/hankel/issues/new?template=feature-request.md>`_.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :-)

Get Started!
------------

Ready to contribute? Here's how to set up ``hankel`` for local development.

1. Fork the ``hankel`` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/hankel.git

3. Install your local copy into a virtualenv. Assuming you have
   virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv hankel
    $ cd hankel/
    $ pip install -e .[dev]

   Note the optional extra install of the development dependencies.

   If you are using ``conda``, setup your environment the usual way, i.e.
   ``conda create -n hankel`` before installing as above.

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox::

    $ flake8 hankel tests
    $ py.test
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

   $ git add .
   $ git commit -m "Your detailed description of your changes."
   $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated.
   Ensure your new functionality is documented in the relevant docstring (we use
   Numpy-style docstrings). If it is a significant feature, it will also
   be preferable to add a demo notebook to the ``docs/demos/`` directory.
3. If you implement an important new feature, consider adding the
   feature to the list in README.rst.
4. The pull request should work for Python 2.7, 3.5, 3.6, and 3.7. Check
   https://travis-ci.org/steven-murray/hankel/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

$ py.test tests.test_hankel


Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run::

$ bumpversion patch # possible: major / minor / patch
$ git push
$ git push --tags

Travis will then deploy to PyPI if tests pass.
