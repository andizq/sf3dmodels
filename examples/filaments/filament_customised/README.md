# Dust continuum and line emission for a filament model

Type the following commands:

```bash
$ python filament_customised.py
$ curl https://home.strw.leidenuniv.nl/~moldata/datafiles/co.dat -o co.dat
$ lime -nSG -p 4 rt-lime.c   #Normal output (-n); sf3dmodels & irregular grid (-SG) and 4 threads (-p 4)
$ python make_moment.py
```