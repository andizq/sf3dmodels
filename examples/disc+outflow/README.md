# Disc + outflow: dust and CO emission in LIME, using H2 (disc) and Hplus (outflow) species

Type the following commands:

```bash
python make_outflow.py
python make_disc.py
python overlap_models.py
curl https://home.strw.leidenuniv.nl/~moldata/datafiles/co.dat -o co.dat
lime -nS -p 8 rt-lime19.c
```