MUSE
Copy *ALL* files from the ESO archive. Include in the same directory the eso-rex calibrations files. 
Needs to be edited the name of the astrometry_wcs.
astrometry_wcs='astrometry_wcs_wfm_2015-09-10.fits'

select_muse_xml_v2.py: 

Reads xml files and search for each file listed there. It creates the OBS Number (i.e. the ODB) and copy all the related files to the Number.
the otput is an archive named  copia_todo.sh.:
chmod +x copia_todo.sh
./copia_todo.sh

ordena_v2.py: edit the astroemtry_wcs file.
It read the headers of all fits and fits.fz files and it creates the corresponding  .sof files to be used with esorex.

creates a final script file named REDUCE_MUSE.sh
chmod +x REDUCE_MUSE.sh
./REDUCE_MUSE.sh


