# Dust continuum and line emission from a broken disc, triggered by a binary system, simulated with the Phantom SPH code 

Download the simulation snapshot [here](https://girder.hub.yt/#user/5da06b5868085e00016c2dee/folder/5f9173da68085e0001d27e7e)
Simulation credit: Facchini et al. 2018

Type the following commands in a terminal:

```bash
python read_write_snap.py #Read and clean snapshot; create Voronoi; write files for Polaris
polaris POLARIS.cmd #Compute dust temperature distribution using two radiation sources
python read_temp.py #Read Polaris temperatures and provide Lime with the final file
curl https://home.strw.leidenuniv.nl/~moldata/datafiles/co.dat -o co.dat #Download molecule info
lime -nS -p 8 rt-model-co.c #Run Lime in sf3dmodels mode to compute line and continuum emission
```