<common>
	<dust_component>	"/path/to/polaris/input/dust/silicate.dat" "plaw" 0.625 3500.0 5e-09 2.5e-07 -3.5
	<dust_component>	"/path/to/polaris/input/dust/graphite_perpend.dat" "plaw" 0.25 2250 5e-09 2.5e-07 -3.5
	<dust_component>	"/path/to/polaris/input/dust/graphite_parallel.dat" "plaw" 0.125 2250 5e-09 2.5e-07 -3.5
	
	<axis1>	1	0	0
	<axis2>	0	1	0

	<write_inp_midplanes>	256
	<write_out_midplanes>	256
	<plot_inp_midplanes>	0
	<plot_out_midplanes>	0
	<midplane_zoom>		1

	<mass_fraction>		0.01
	<mu>			2.0

	<phase_function>	PH_MIE

	<nr_threads>		-1
</common>

<task> 1
	<source_star nr_photons = "1000000"> 2.34494430e+11  2.92622314e+11 -6.28616012e+09	2.5	9500
	<source_star nr_photons = "1000000"> -2.32978652e+11 -2.87857614e+11  8.84130821e+09 	2.5	9500

	<cmd>			CMD_TEMP

	<path_grid>		"./voronoi_grid.dat"
	<path_out>		"./"

	<dust_offset>		0

	<adj_tgas>		1
	<sub_dust>		1

	<conv_dens>		1.0
	<conv_len>		1.0
	<conv_mag>		1.0
	<conv_vel>		1.0
</task>

