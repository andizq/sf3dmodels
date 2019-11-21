# Dust (silicate) continuum emission from powerlaw-shells model

Type the following commands:

```bash
python make_shellsplaw.py
radmc3d mctherm
radmc3d image lambda 1000
python plot_png_image.py
```