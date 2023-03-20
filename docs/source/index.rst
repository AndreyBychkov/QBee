.. QBee documentation master file, created by
   sphinx-quickstart on Wed Oct 13 01:22:41 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to QBee's documentation!
================================

**qbee** is a package for transforming differential equations into polynomial and quadratic polynomial forms via introducing new variables.
It includes functionality for:

* exact transformation of Ordinary Differential Equations (ODEs) into a form with a polynomial right-hand side (**polynomialization**).
* exact transformation of polynomial ODEs intro a form with right-hand side being quadratic polynomials (**quadratization**).
* obtaining a general (vectorized) form of new variables to quadratize linearly coupled polynomial ODEs given in vectorized form (**dimension-agnostic quadratization**)
* automatized transformation of polynomialized and quadratized systems of ODEs into a form suitable for **numerical integrators**.

Background - Why we created **qbee**
------------------------------------

Lorem ipsum add articles and examples


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   quadratization
   inputs_and_parameters
   polynomialization
   dimension_agnostic
   numerical_integration
   articles_and_resources
   api