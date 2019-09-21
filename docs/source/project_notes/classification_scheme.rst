Classification Scheme
=====================


Optimizing the Figure of Merit (FOM)
------------------------------------

Targets are classified based on their position in the following phase
space:

.. math::

    x \def \chi^2_{blue}(Ia) - \chi^2_{blue}(91bg)
    y \def \chi^2_{red}(Ia) - \chi^2_{red}(91bg)

To determine where in this phase space we should draw the boundaries
separating each classification, we define the following values for a given type
of object (e.g. 91bg-like objects):

 - :math:`N_{total}` is the total number of objects with the given type (i.e. the "truth")
 - :math:`N_{true}` is the number object correctly classified as the given type
 - :math:`N_{false}` is the number of objects falsely classified the given type

We then maximize the FOM:

.. math::

    FOM = \frac{N_{true}}{N_{total}} \times \frac{N_{true}}{N_{true} + N_{false}}


