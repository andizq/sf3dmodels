# Dust (silicate) continuum emission from powerlaw-shells model

Type the following commands:

```bash
python make_shellsplaw.py #Writes the input files for radmc3d
cp path/to/radmc-3d/version_0.41/examples/run_1dpp_dust/dustkappa_silicate.inp ./ #Copies file with absorption properties of silicates into folder
radmc3d mctherm #Computes dust temperature distribution using thermal MonteCarlo simulation
radmc3d image lambda 1000 #Dust continuum image
python plot_png_image.py  #Plotting output
```