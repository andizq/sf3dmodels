# Dust continuum and line emission from a composite model: filament + core

Type the following commands:

```bash
$ python filament_default.py
$ python shellsplaw_infall.py
$ python overlap_submodels.py
$ curl https://home.strw.leidenuniv.nl/~moldata/datafiles/co.dat -o co.dat
$ lime -nS -p 4 rt-lime.c   #Normal output (-n); sf3dmodels & regular grid (-S) and 4 threads (-p 4)
$ python make_moment.py
```