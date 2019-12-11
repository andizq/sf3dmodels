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
	<detector_dust nr_pixel = "256*256">	0.00085	0.00085	1	1	0.0	0.0	3.086e+19

	<source_background nr_photons = "100000.0">	0.0	0.0	0.0	0.0	0.0	0.0	0.0

	<cmd>			CMD_DUST_EMISSION
	<path_grid>		"./grid_temp.dat"
	<path_out>		"./polaris_dustcont/"

	<max_subpixel_lvl>	1
	<f_highJ>		0.25
	<f_c>			0.0

	<conv_dens>		1.0
	<conv_len>		1.0
	<conv_mag>		1.0
	<conv_vel>		1.0
</task>

