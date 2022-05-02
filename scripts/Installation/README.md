# Quick tips for using PSP

1. On Gaanam clusters, you can use the python version at:
```/home/modules/anaconda3/bin/python3```

2. Add these to your ```~/.bashrc```

```
# PACKMOL PATH
export PACKMOL_EXEC='/home/appls/utility/packmol/packmol'

# pysimm PATHs
export PYTHONPATH=$PYTHONPATH:'/home/appls/utility/pysimm'
export PATH=$PATH:'/home/appls/utility/pysimm/bin'

# AmberTools PATH
export ANTECHAMBER_EXEC='/home/modules/anaconda3/bin/antechamber'
```