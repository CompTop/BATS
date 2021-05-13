.. installation instructions:

Installation
============

You can obtain BATS from GitHub:

.. code-block:: console

   $ git clone https://github.com/bnels/BATS.git

Dependencies
------------

BATS is written in C++ and requires a compiler which is able to
understand C++17 syntax.  It has been run and tested on Linux
(Ubuntu and Fedora distributions) using the gcc compiler.

BATS itself is header-only, and is built using STL data structures,
meaning you don't need to find other headers.

To get the most out of BATS, you should compile with OpenMP.

.. code-block:: console

   $ g++ ...<other flags> -std=c++17 -fopenmp


Testing Installation
--------------------

Once you have BATS on your computer, you can test your installation.
from the root of your repository:

.. code-block:: console

   $ cd tests
   $ make test -j
