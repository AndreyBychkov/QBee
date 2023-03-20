API Reference
=============

This section contains a summary of QBee functions and variables.

.. automodule:: qbee
   :members: polynomialize_and_quadratize, quadratize, polynomialize, quadratize_dimension_agnostic,
               pruning_by_nodes_processed, pruning_by_vars_number, pruning_by_declining_variables,
               pruning_by_elapsed_time, pruning_by_nodes_without_quadratization_found, default_pruning_rules, without_variables,
               functions, parameters,
               print_qbee, str_qbee,
               default_generation, generation_semidiscretized, default_scoring, smd_scoring, aeqd_scoring,
   :member-order: bysource

.. data:: qbee.INDEPENDENT_VARIABLE

   Independent variable for functions.

   Example:
      >>> from qbee import *
      >>> from qbee.experimental import to_odeint
      >>> from sympy import sin
      >>> x, u = functions("x, u")
      >>> res = polynomialize_and_quadratize([(x, sin(x) + u)])
      >>> t = INDEPENDENT_VARIABLE
      >>> my_odeint = to_odeint(res, {x: 1}, inputs={u: sin(t)})