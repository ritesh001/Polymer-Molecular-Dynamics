---
sidebar_label: ForceField
title: core.ForceField
---

## OPLS Objects

```python
class OPLS(ForceField)
```

Template OPLS object to contain OPLS force field settings

**Attributes**:

- `charge_method` _str_ - Charge method; has to be one of `"cm1a-lbcc"` or
  `"cm1a"`; default: `"cm1a-lbcc"`

## GAFF2 Objects

```python
class GAFF2(ForceField)
```

Template GAFF2 object to contain GAFF2 force field settings

**Attributes**:

- `charge_method` _str_ - Charge method; has to be one of `"gasteiger"` or
  `"am1bcc"`; default: `"gasteiger"`
