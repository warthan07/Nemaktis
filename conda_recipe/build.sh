##############
# BPM solver #
##############

mkdir -p bpm/bin
cd bpm/bin
rm -rf *

cmake ..
make install python-install
cd ../..


##############
# RTM solver #
##############

if [ -d rtm ]; then
	mkdir -p rtm/bin
	cd rtm/bin
	rm -rf *
	
	cmake ..
	make install
	cd ../..
fi


###############
# DTMM solver #
###############

cd dtmm
python setup.py install
cd ..


########################
# High-level interface #
########################

cd nemaktis
python setup.py install
