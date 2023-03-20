Installation
=================

QBee can be installed on virtually any computer with Python.
Python can be installed via `official website <https://www.python.org/downloads/>`_
and QBee officially supports Python 3.9 and above.

Pip / PyPI
-----------------

Official distribution of QBee is located on `PyPI <https://pypi.org/project/qbee/>`_
and can be installed via pip in your console:

.. code-block:: console

   $ pip install qbee

Git
-------------

If you want to contribute to QBee or like to get each possible version of QBee, install it from git.
To download repository, execute the following from the command line:

.. code-block:: console

    $ git clone https://github.com/AndreyBychkov/QBee.git --branch master

You can pick a different branch if you need it.
For instance, ``dev`` branch contains latest updates, and others and contain unreleased new features.

To update to the latest version, go into your repository and execute:

.. code-block:: console

    $ git pull origin master

If you want to install QBee, you can do it from the repository by using pip:

.. code-block:: console

    $ pip install .


Run QBee
------------------

After installation, it is best to verify that your the freshly-installed package works.
To do this, start up Python and import the QBee library:

.. code-block:: console

    $ python
    >>> from qbee import *

From here, calculate a simple quadratization:

.. code-block:: python

    >>> x = functions("x")
    >>> quadratize([(x, x**5)]).print()
    Introduced variables:
    w0 = x**4

    x' = w0*x
    w0' = 4*w0**2
