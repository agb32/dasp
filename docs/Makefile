# Makefile for aosim documentation

# \/ start of user modifiable area

PYTHON = python
   # \/ add your document filename stem here to get it made
DOCUMENTS = introduction dummies spatial_filter overview paramgui daspctrl schema daspsetup simnews conjugation    setupguirequirements simapi-orig
#fdpcg cellprocessor fpgaMVM wfsPipe

# /\ end of user modifiable area

all:  ${DOCUMENTS:=.pdf}
	echo ">>>"$<
	${PYTHON} makeLatex.py make
	(cd modules/ && latex api.tex && latex api.tex && makeindex api && latex api.tex && latex api.tex && dvipdf api)
	cp modules/api.pdf .
	#echo "Available documentation in PDF format (use acroread):"
	#ls *.pdf

.INTERMEDIATE: ${DOCUMENTS:=.bbl} ${DOCUMENTS:=.log} ${DOCUMENTS:=.aux} ${DOCUMENTS:=.ilg} ${DOCUMENTS:=.idx} ${DOCUMENTS:=.dia} ${DOCUMENTS:=.out} ${DOCUMENTS:=.toc} ${DOCUMENTS:=.ind}

%.dvi %.log %.ind %.idx %.ilg %.aux %.dia %.out %.toc %.bbl: %.tex
	latex $<
	latex $<
	makeindex $*
	if [ -e $*.bib ] ; then bibtex $* ; fi
	latex $<
	latex $<

%.ps: %.dvi
	dvips -Ppdf -o $@ $<

%.pdf: %.dvi
	dvipdf $<

clean: 
	rm -rf modules
	rm -f ${DOCUMENTS:=.pdf} simapi.tex *.*~

ubuntu1404:
	sudo apt-get install texlive-latex-base texlive-latex-recommended python-epydoc



# \/ below are older, explicit rules

#all:  ${DOCUMENTS:=.pdf} dummiesGuide.pdf overview.pdf simctrl.pdf paramgui.pdf schema.pdf setupGUIrequirements.pdf wfsPipe.pdf fpgaMVM.pdf simsetup.pdf simnews.pdf conjugation.pdf fdpcg.pdf cellprocessor.pdf
#	${PYTHON} makeLatex.py make
#	echo "Available documentation in PDF format (use acroread):"
#	ls *.pdf
#
#tomographicSimulations.pdf: tomographicSimulations.tex
#	latex tomographicSimulations.tex
#	latex tomographicSimulations.tex
#	makeindex tomographicSimulations
#	bibtex tomographicSimulations
#	latex tomographicSimulations.tex
#	latex tomographicSimulations.tex
#	dvips -Ppdf -o tomographicSimulations.ps tomographicSimulations.dvi
#	ps2pdf13 tomographicSimulations.ps
#
#cellprocessor.pdf: cellprocessor.tex
#	latex cellprocessor.tex
#	latex cellprocessor.tex
#	makeindex cellprocessor
#	bibtex cellprocessor
#	latex cellprocessor.tex
#	latex cellprocessor.tex
#	dvips -Ppdf -o cellprocessor.ps cellprocessor.dvi
#	ps2pdf13 cellprocessor.ps
#
#conjugation.pdf: conjugation.tex
#	latex conjugation.tex
#	latex conjugation.tex
#	makeindex conjugation
#	latex conjugation.tex
#	latex conjugation.tex
#	dvips -Ppdf -o conjugation.ps conjugation.dvi
#	ps2pdf13 conjugation.ps
#
#simnews.pdf: simnews.tex
#	latex simnews.tex
#	latex simnews.tex
#	makeindex simnews
#	bibtex simnews
#	latex simnews.tex
#	latex simnews.tex
#	dvips -Ppdf -o simnews.ps simnews.dvi
#	ps2pdf13 simnews.ps
#
#fdpcg.pdf: fdpcg.tex
#	latex fdpcg.tex
#	latex fdpcg.tex
#	makeindex fdpcg
#	bibtex fdpcg
#	latex fdpcg.tex
#	latex fdpcg.tex
#	dvips -Ppdf -o fdpcg.ps fdpcg.dvi
#	ps2pdf13 fdpcg.ps
#
#simsetup.pdf: simsetup.tex
#	latex simsetup.tex
#	latex simsetup.tex
#	makeindex simsetup
#	bibtex simsetup
#	latex simsetup.tex
#	latex simsetup.tex
#	dvips -Ppdf -o simsetup.ps simsetup.dvi
#	ps2pdf13 simsetup.ps
#
#fpgaMVM.pdf: fpgaMVM.tex
#	latex fpgaMVM.tex
#	latex fpgaMVM.tex
#	makeindex fpgaMVM
#	bibtex fpgaMVM
#	latex fpgaMVM.tex
#	latex fpgaMVM.tex
#	dvips -Ppdf -o fpgaMVM.ps fpgaMVM.dvi
#	ps2pdf13 fpgaMVM.ps
#
#wfsPipe.pdf: wfsPipe.tex
#	latex wfsPipe.tex
#	latex wfsPipe.tex
#	makeindex wfsPipe
#	bibtex wfsPipe
#	latex wfsPipe.tex
#	latex wfsPipe.tex
#	dvips -Ppdf -o wfsPipe.ps wfsPipe.dvi
#	ps2pdf13 wfsPipe.ps
#
#schema.pdf: schema.tex
#	latex schema.tex
#	latex schema.tex
#	makeindex schema
#	#bibtex schema
#	latex schema.tex
#	latex schema.tex
#	dvips -Ppdf -o schema.ps schema.dvi
#	ps2pdf13 schema.ps
#
#simapi.pdf:simapi.tex
#	latex simapi.tex
#	latex simapi.tex
#	makeindex simapi
#	bibtex simapi
#	latex simapi.tex
#	latex simapi.tex
#	dvips -Ppdf -o simapi.ps simapi.dvi
#	ps2pdf13 simapi.ps
#
#myeso: useeso.tex
#	latex useeso.tex
#	latex useeso.tex
#	dvips -Ppdf -o useeso.ps useeso
#	ps2pdf useeso.ps
#skeleton.ps: skeleton.tex
#	latex skeleton.tex
#	dvips -Ppdf -o skeleton.ps skeleton
#	ps2pdf skeleton.ps
#eso: esotemplate.tex
#	latex esotemplate.tex
#	latex esotemplate.tex
#	dvips -Ppdf -o eso.ps esotemplate
#	ps2pdf eso.ps
#overview.pdf: overview.tex
#	latex overview.tex
#	latex overview.tex
#	makeindex overview
#	bibtex overview
#	latex overview.tex
#	latex overview.tex
#	dvips -Ppdf -o overview.ps overview.dvi
#	ps2pdf13 overview.ps
#
#dummiesGuide.pdf: dummies.tex
#	latex dummies.tex
#	latex dummies.tex
#	makeindex dummies
#	bibtex dummies
#	latex dummies.tex
#	latex dummies.tex
#	dvips -Ppdf -o dummiesGuide.ps dummies.dvi
#	ps2pdf13 dummiesGuide.ps
#paramgui.pdf: paramgui.tex
#	latex paramgui.tex
#	latex paramgui.tex
#	makeindex paramgui
#	bibtex paramgui
#	latex paramgui.tex
#	latex paramgui.tex
#	dvips -Ppdf -o paramgui.ps paramgui.dvi
#	ps2pdf13 paramgui.ps
#simctrl.pdf: simctrl.tex
#	latex simctrl.tex
#	latex simctrl.tex
#	makeindex simctrl
#	bibtex simctrl
#	latex simctrl.tex
#	latex simctrl.tex
#	dvips -Ppdf -o simctrl.ps simctrl.dvi
#	ps2pdf13 simctrl.ps
#setupGUIrequirements.pdf: setupguirequirements.tex
#	latex setupguirequirements.tex
#	latex setupguirequirements.tex
#	makeindex setupguirequirements
#	bibtex setupguirequirements
#	latex setupguirequirements.tex
#	latex setupguirequirements.tex
#	dvips -Ppdf -o setupGUIrequirements.ps setupguirequirements.dvi
#	ps2pdf13 setupGUIrequirements.ps
