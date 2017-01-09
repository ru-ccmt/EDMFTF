include Makefile.in

SRC = impurity/ctqmc impurity/atom_ctqmc impurity/atomg impurity/broad impurity/NCA impurity/OCA impurity/atomd impurity/fort impurity/maxent_source impurity/ancont impurity/skrams optics magnetSusc/dmft_bubble magnetSusc/dmftgk downfold putils/MMoment putils/local_axes putils/wavef exactDC dmft1 dmft2 dftKS dftso


all:
	for dir in ${SRC} - ; do \
		if [ "$$dir" = "-" ] ; then break ; fi ; \
		(cd $$dir && ${MAKE} all ; ) ; done

clean:
	for dir in ${SRC} - ; do \
		if [ "$$dir" = "-" ] ; then break ; fi ; \
		(cd $$dir && ${MAKE} clean ; ) ; done
	rm -f *.pyc


install :
	if [ ! -d ${DESTDIR} ] ; then mkdir -p ${DESTDIR} ; fi
	if [ ! -d ${DESTDIR}/database ] ; then mkdir -p ${DESTDIR}/database ; fi
	cp dmft1/dmft ${DESTDIR}/
	cp dmft1/*.so ${DESTDIR}/
	cp dmft1/chemicalP.py ${DESTDIR}/
	cp impurity/fort/*.so ${DESTDIR}/
	cp impurity/database/* ${DESTDIR}/database/
	cp impurity/NCA/nca ${DESTDIR}/
	if [ ${wildcard impurity/OCA/oca} ] ; then cp impurity/OCA/oca ${DESTDIR}/; fi
	cp impurity/atomg/atom ${DESTDIR}/
	cp impurity/broad/broad ${DESTDIR}/
	cp impurity/atomd/*.so ${DESTDIR}/
	cp impurity/atomd/atom_d.py ${DESTDIR}/
	cp impurity/atom_ctqmc/atomh ${DESTDIR}/
	cp impurity/ctqmc/ctqmc ${DESTDIR}/
	cp impurity/*.py ${DESTDIR}/
	cp impurity/link/link.py ${DESTDIR}/
	cp dmft2/dmft2 ${DESTDIR}/
	cp dmft2/fermif ${DESTDIR}/
	cp dmft2/combineud ${DESTDIR}/
	cp python/*.py ${DESTDIR}/
	cp optics/dmftopt ${DESTDIR}/
	cp downfold/*.py ${DESTDIR}/
	cp downfold/*.so ${DESTDIR}/
	cp impurity/maxent_source/maxent_routines.so ${DESTDIR}/
	cp impurity/maxent_source/maxentropy.py ${DESTDIR}/
	cp impurity/maxent_source/maxent_run.py ${DESTDIR}/
	cp impurity/skrams/skrams ${DESTDIR}/
	cp magnetSusc/dmft_bubble/dmft_real_bubble ${DESTDIR}/
	cp magnetSusc/dmftgk/dmftgk ${DESTDIR}/
	cp magnetSusc/python/*.py ${DESTDIR}/
	cp putils/findDiagonal/*.py ${DESTDIR}/
	cp putils/MMoment/*.so ${DESTDIR}/
	cp putils/MMoment/local2global.py ${DESTDIR}/
	cp putils/MMoment/MagneticMoment.py ${DESTDIR}/
	cp putils/MMoment/Wigner.py ${DESTDIR}/
	cp putils/MMoment/transformations.py ${DESTDIR}/
	cp putils/local_axes/localaxes.py ${DESTDIR}/
	cp putils/local_axes/w2k_nn.so ${DESTDIR}
	cp putils/wavef/w2k_atpar.so ${DESTDIR}
	cp putils/wavef/wavef.py ${DESTDIR}
	cp dftKS/lapw1{c,} ${DESTDIR}
	cp dftso/lapwso ${DESTDIR}
	cp exactDC/*.py ${DESTDIR}
	cp exactDC/*.so ${DESTDIR}

