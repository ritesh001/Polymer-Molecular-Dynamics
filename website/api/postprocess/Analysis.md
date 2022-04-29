---
sidebar_label: Analysis
title: postprocess.Analysis
---

#### calculate\_Tg

```python
def calculate_Tg(result_fname: str, make_plot: bool = True) -> int
```

Method to calculate glass transition temperature based on the
result file obtained from TgMeasurement Procedure

**Arguments**:

- `result_fname` _str_ - Name of the result file from TgMeasurement
  Procedure
- `make_plot` _bool_ - Whether to make a plot to visualize the fitting
  

**Returns**:

- `Tg` _int_ - Glass transition temperature of the system

