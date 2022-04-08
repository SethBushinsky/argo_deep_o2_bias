# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# ## Functions fo carbonate system calculations (pH, TALK etc)

import numpy as np
import glob, os
import pandas as pd
import xarray as xr
import gsw
import PyCO2SYS as pyco2
