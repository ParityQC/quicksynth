<!-- Copyright (c) 2026 Parity Quantum Computing GmbH
     SPDX-License-Identifier: MIT -->

# API reference

```{eval-rst}
.. automodule:: quicksynth
   :no-members:
```

## Synthesis

```{eval-rst}
.. currentmodule:: quicksynth

.. autofunction:: synthesise

.. autofunction:: synthesise_rotation

.. autoclass:: QuickSynth

.. autoclass:: SynthesisResult
```

## Settings

```{eval-rst}
.. currentmodule:: quicksynth

.. autoclass:: SynthesisConfig
   :exclude-members: __init__, __new__
```

## Targets

```{eval-rst}
.. currentmodule:: quicksynth

.. autofunction:: rotation

.. autofunction:: nearest_unitary
```

## Exact layer

```{eval-rst}
.. automodule:: quicksynth.api
   :no-members:

.. currentmodule:: quicksynth

.. autoclass:: CyclotomicMatrix

.. autofunction:: word_unitary

.. autofunction:: price

.. autofunction:: matrix_word

.. autofunction:: normal_form

.. autofunction:: bloch_image
```
