<common>
	<dust_component>	"/Users/andrespipecar42/polaris/input2/dust/silicate.dat" "plaw" 0.625 3500.0 5e-09 2.5e-07 -3.5
	<dust_component>	"/Users/andrespipecar42/polaris/input2/dust/graphite_perpend.dat" "plaw" 0.25 2250 5e-09 2.5e-07 -3.5
	<dust_component>	"/Users/andrespipecar42/polaris/input2/dust/graphite_parallel.dat" "plaw" 0.125 2250 5e-09 2.5e-07 -3.5
	
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
	<source_star nr_photons = "10000000"> 0  0  0	0.1  3766

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

