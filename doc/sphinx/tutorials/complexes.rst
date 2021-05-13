Complexes
=========

This tutorial will cover building and using complexes which represent
topological spaces.

Interfaces
----------

Chain Functor
^^^^^^^^^^^^^
There are several functions a complex should implement in order to
be compatible with the chain functor.

Filtrations
^^^^^^^^^^^
Filtrations can be created from a complex as well as a list of filtration
values.  However, in order to construct filtrations incrementally the
following methods should be implemented.
