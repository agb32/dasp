#!/usr/bin/env python
import time
import os
import string
import util.parseSimXml
def makesim():
    print("""Welcome to daspbuilder.  This will help you to set up simple simulations.
This serves as an entry point into the functionality of dasp.  If you require
more complex simulations, you are able to add these using the daspsetup.py tool
and by editing parameter files.  For example, you may wish to add modal or
tip-tilt DMs, or you may wish to move some DMs to be open-loop
and others closed-loop.  All this functionality is possible if you dig deep
enough.

""")
    simtyp=raw_input("""Please enter your simulation type:
[1] SCAO (including XAO)
 2  GLAO
 3  LTAO
 4  MCAO
 5  MOAO
 6  Interaction matrix
 7  Wavefront recorder
 8  MPI instructions
 9  Help and notes
""")
    if len(simtyp)==0:
        simtyp="1"
    simtyp=int(simtyp)
    if simtyp==1:
        makescao()
    elif simtyp==2:
        makemcao(ndm=1,defaultdirname="glaoSim",fname="glao")
    elif simtyp==3:
        makeltao()
    elif simtyp==4:
        makemcao()
    elif simtyp==5:
        makemoao()
    elif simtyp==6:
        makepmx()
    elif simtyp==7:
        makelearn()
    elif simtyp==8:
        printmpi()
    elif simtyp==9:
        printhelp()
    else:
        print("Idiot!")

def printhelp():
    print("""General help and guidance:

This tool will aid you in setting up simple AO simulations.  After
this tool has been used, you should then edit the params.py file, to
match the simulation to your requirements, for example, telescope
diameter, number of sub-apertures, atmospheric turbulence
characteristics, etc.  However, please note that if you change the
number of DMs in the parameter file, in most cases, you will need to
use daspbuilder.py to regenerate a corresponding simulation (or add
the DMs yourself using daspsetup.py).

The performance of the wide-field AO systems simulated here is
generally fairly poor when using the parameter file as it is:
Correction on a 4.2m telescope using 70cm sub-apertures, with LGS on a
2 arcminute diameter ring, and NGS on a 3 arcminute diameter ring, is
never going to be good.  If you want to convince yourself that the
simulation is working, I suggest reducing the asterism diameters.
Alternatively, daspsetup.py can be used to add an uncorrected science
source, which can then be used for comparison.

To view the simulation configurations, use:
daspsetup.py SIMFILE.xml  (where SIMFILE is the generated filename)

To connect to a running simulation (to view PSFs, WFSs, DMs, etc in action) use:
daspctrl.py

To get current Strehl, etc, from the commandline, use
daspanalyse.py
""")

def printmpi():
    print("""Each simulation gives you the option of generating an MPI version.
In general, to run this, instead of using python, you need to call:
mpirun -np XXX python FILE.py PARAMS.py
where XXX is the number of MPI processes (given in the top line of FILE.py),
FILE.py is your simulation file, and PARAMS.py is the parameter file.

All the standard MPI options are available, e.g. -machinefile, -hostlist, etc
(though depending on your MPI installation, this may vary).  You may also need
to give a full path to python, or even mpipython, depending on which Python MPI
module you are using.
""")

def save(txt,dirname,name,doexit=1):
    fname=os.path.join(dirname,name)
    if os.path.exists(fname):
        yn=raw_input("Overwrite existing file %s? [y] "%fname)
        if len(yn)>0 and yn!="y":
            if doexit:
                print "Exiting"
                sys.exit(0)
            else:
                print "Not saving %s"%fname
        else:#make a backup.
            os.rename(fname,fname+".bak")
    open(fname,"w").write(txt)

def makescao():
    nsci=raw_input("Please enter the number of science targets required\n(if you choose 1, it will be on-axis)\n[1]: ")
    if len(nsci)==0:
        nsci=1
    else:
        nsci=int(nsci)
    if nsci==1:
        xmltxt=scao1sciTxt
        p=util.parseSimXml.parseSimXml(xmltxt)
        pystr=p.makePython()
    else:
        gr=""
        for i in range(nsci):
            gr+="sci%d,"%(i+1)
        xmltxt=scaoNsciTxt.replace("sci1",gr)
        p=util.parseSimXml.parseSimXml(xmltxt)
        pystr=p.makePython()
    dirname="scaoSim"#%time.strftime("%y%m%d_%H%M%S")
    tmp=raw_input("\n\nPlease enter an output directory name:\n[%s] "%dirname)
    if len(tmp)>0:
        dirname=tmp
    if os.path.exists(dirname):
        print("Using existing directory %s"%dirname)
    else:
        os.mkdir(dirname)
    pyname="scao.py"
    xmlname="scao.xml"
    paramname="params.py"
    save(pystr,dirname,pyname)
    save(xmltxt,dirname,xmlname)
    save(scaoParamsTxt.replace("nsci=1","nsci=%d"%nsci),dirname,paramname)
    print("""Simulation generated.
To run in %s/, use:
python scao.py params.py
To view the simulation configuration, use:
daspsetup.py scao.xml
To connect to a running simulation (to view PSFs, WFSs, DMs, etc in action) use:
daspctrl.py
or:
daspanalyse.py

Please edit params.py to alter the simulation to meet your requirements.

Some notes about this simulation:
It will first of all take an interaction matrix.
This will then be inverted to give a control matrix, using a simple SVD.
This least-squares wavefront reconstructor will then be used for AO correction.

Several improvements are possible (left as an excercise for the user):
Minimum variance wavefront reconstruction from the interaction matrix.
Model-based reconstruction (e.g. CuReD, SOR, etc)

If you already have a control matrix, and thus do not need to generate one, you
can use --user=nopoke commandline option

"""%dirname)



def makemcao(nlgs=None,nngs=None,nsci=None,ndm=None,mpi=None,dirname=None,defaultdirname="mcaoSim",fname="mcao",dopoke=1,addvdm=0):
    if nlgs==None:
        nlgs=raw_input("Please enter the number of LGS required: [4] ")
        if len(nlgs)==0:
            nlgs=4
        else:
            nlgs=int(nlgs)
    if nngs==None:
        nngs=raw_input("Please enter the number of NGS required: [4] ")
        if len(nngs)==0:
            nngs=4
        else:
            nngs=int(nngs)
    if nsci==None:
        nsci=raw_input("Please enter the number of science PSFs required: [1] ")
        if len(nsci)==0:
            nsci=1
        else:
            nsci=int(nsci)
    if ndm==None:
        ndm=raw_input("Please enter the number of DMs required: [3] ")
        if len(ndm)==0:
            ndm=3
        else:
            ndm=int(ndm)
    if mpi==None:
        mpi=raw_input("Use MPI?: y/[n] ")
        if len(mpi)==0 or mpi!="y":
            mpi=0
        else:
            mpi=1
    if dirname==None:
        tmp=raw_input("\n\nPlease enter an output directory name:\n[%s] "%defaultdirname)
        if len(tmp)>0:
            dirname=tmp
        else:
            dirname=defaultdirname
        if os.path.exists(dirname):
            print("Using existing directory %s"%dirname)
        else:
            os.mkdir(dirname)


    ngs=nngs+nlgs
    txt=""
    for i in range(ngs):
        txt+="%d,"%(i+1)
    if mpi:
        xmltxt=mcaoTxt.replace('groupShare idstr="1" cpu="[]"','groupShare idstr="%s" cpu="[%s]"'%(txt,txt))
    else:#no MPI
        xmltxt=mcaoTxt.replace('groupShare idstr="1"','groupShare idstr="%s"'%txt)
    txt=""
    mpitxt=""
    for i in range(nsci):
        txt+="sci%d,"%(i+1)
        mpitxt+="%d,"%(i+1+ngs)
    if mpi:
        xmltxt=xmltxt.replace('groupShare idstr="sci1" cpu="[]"','groupShare idstr="%s" cpu="[%s]"'%(txt,mpitxt))
    else:
        xmltxt=xmltxt.replace('groupShare idstr="sci1"','groupShare idstr="%s"'%txt)

    if dopoke==0:
        xmltxt=xmltxt.replace("ctrl.doInitialPokeThenRun()","pass")


    #now add the DMs.  Need to:
    #Move tomoRecon, wfscent, science and grouping down.
    #Remove existing DM
    #Break connections
    #Add new DMs
    objlist=xmltxt.split("</simulationObject>")
    newlist=[]
    olddmlist=[]
    gsdmlist=[]
    scidmlist=[]
    for obj in objlist:
        if "304" in obj:#y position of wfscent and science
            obj=obj.replace("304","%d"%(224+80*ndm))
        elif "384" in obj:
            obj=obj.replace("384","%d"%(304+80*ndm))
        if "<groupShare" in obj:
            obj=obj.replace("344","%d"%(264+80*ndm))
        if "tomoRecon" in obj:
            #adjust the connectto...
            txt=""
            for i in range(ndm):
                txt+="%d, "%(i+16)
            for i in range(ndm):
                txt+="%d, "%(i+16+ndm)
            obj=obj.replace('connectto="[5, 13]"','connectto="[%s]"'%txt)
            txt=""
            for i in range(ndm):
                txt+="[(129, %d), (129, %d)],\n"%(417+(ndm-1)*80,192+i*80)
            obj=obj.replace("[(129, 417), (129, 192)],\n",txt)
        if "iatmos" in obj:
            if 'tag="3"' in obj:
                obj=obj.replace('connectto="[5]"','connectto="[16]"')
            elif 'tag="11"' in obj:
                obj=obj.replace('connectto="[13]"','connectto="[%d]"'%(16+ndm))
        if "wfscent" in obj:
            obj=obj.replace('connectfrom="[5]"','connectfrom="[%d]"'%(16+ndm-1))
        if "science.science" in obj:
            obj=obj.replace('connectfrom="[13]"','connectfrom="[%d]"'%(16+2*ndm-1))

        if "science.xinterp_dm" not in obj:
            newlist.append(obj)
        else:
            olddmlist.append(obj)
    dm=olddmlist[0]
    for i in range(ndm):
        txt=dm.replace("80,224","80,%d"%(224+i*80))
        txt=txt.replace('tag="5"','tag="%d"'%(16+i))
        if i==ndm-1:
            connto=7
        else:
            connto=16+i+1
        txt=txt.replace('connectto="[7]"','connectto="[%d]"'%(connto))
        if i==0:
            connfrom=3
        else:
            connfrom=16+i-1
        txt=txt.replace('connectfrom="[3','connectfrom="[%d'%(connfrom))
        txt=txt.replace("[(129, 417), (129, 192)]","[(129, %d), (129, %d)]"%(417+(ndm-1)*80,192+i*80))
        txt=txt.replace("<sharedTo>\n[13]","<sharedTo>\n[%d]"%(16+ndm+i))
        txt=txt.replace("dm$","dm%dpath$"%i)
        gsdmlist.append(txt)

        txt=dm.replace("80,224","180,%d"%(224+i*80))
        txt=txt.replace('tag="5"','tag="%d"'%(16+i+ndm))
        if i==ndm-1:
            connto=15
        else:
            connto=16+ndm+i+1
        txt=txt.replace('connectto="[7]"','connectto="[%d]"'%(connto))
        if i==0:
            connfrom=11
        else:
            connfrom=16+ndm+i-1
        txt=txt.replace('connectfrom="[3','connectfrom="[%d'%(connfrom))
        txt=txt.replace("[(129, 417), (129, 192)]","[(129, %d), (129, %d)]"%(417+(ndm-1)*80,192+i*80))
        txt=txt.replace("<sharedTo>\n[13]","<sharedTo>\n[]")
        txt=txt.replace("dm$","dm%dpath$"%i)
        scidmlist.append(txt)
        
    objlist=newlist[:2]+gsdmlist+newlist[2:5]+scidmlist+newlist[5:]
    xmltxt=string.join(objlist,"</simulationObject>")
    save(xmltxt,dirname,fname+".xml")
    p=util.parseSimXml.parseSimXml(xmltxt)
    pystr=p.makePython()
    save(pystr,dirname,fname+".py")
    txt=mcaoParamsTxt.replace("nsci=1","nsci=%d"%nsci)
    txt=txt.replace("nlgs=4","nlgs=%d"%nlgs)
    txt=txt.replace("nngs=3","nngs=%d"%nngs)
    txt=txt.replace("ndm=3","ndm=%d"%ndm)
    if addvdm:#for ltao
        txt=txt.replace("dmOverview=dmOverview(dmInfoList,atmosGeom)","""dmInfoList.append(dmInfo('vdm',["sci1"],0,nAct,minarea=0.1,actuatorsFrom=["dm%dpath"%x for x in range(ndm)],maxActDist=1.5,decayFactor=0.95,reconLam=lgsLam))
dmOverview=dmOverview(dmInfoList,atmosGeom)""")
    save(txt,dirname,"params.py")
    if mpi:
        execstr="mpirun -np %d python %s.py params.py"%(ngs+nsci,fname)
        mpitxt="\nThe order of module execution may need editing to maximise performance:\nPlease edit lines starting with 'execOrder' in %s.py to do this.\n(sorry: the logic is difficult to get right for all situations).\n"%fname
    else:
        execstr="python %s.py params.py"%fname
        mpitxt=""

    print("""Simulation generated.
To run in %s/, use:
%s
To view the simulation configuration, use:
daspsetup.py %s.xml
To connect to a running simulation (to view PSFs, WFSs, DMs, etc in action) use:
daspctrl.py
or:
daspanalyse.py

Please edit params.py to alter the simulation to meet your requirements.

Some notes about this simulation:
It will first of all take an interaction matrix.
This will then be inverted to give a control matrix, using a simple SVD.
This least-squares wavefront reconstructor will then be used for AO correction.

The LGS tip-tilt signal is assumed valid, i.e. is used.  

Several improvements are possible (left as an excercise for the user):
Minimum variance wavefront reconstruction from the interaction matrix.
Ignore the LGS tip-tilt signal, in the user-generated control matrix.
Model-based reconstruction (e.g. FrIM, FEWHA, etc)

If you already have a control matrix, and thus do not need to generate one, you
can use --user=nopoke commandline option
%s"""%(dirname,execstr,fname,mpitxt))





def makepmx(nlgs=None,nngs=None,ndm=None,mpi=None,dirname="pokeSim",writeParams=1,fname="poke"):
    if nlgs==None:
        nlgs=raw_input("Please enter the number of LGS required: [4] ")
        if len(nlgs)==0:
            nlgs=4
        else:
            nlgs=int(nlgs)
    if nngs==None:
        nngs=raw_input("Please enter the number of NGS required: [4] ")
        if len(nngs)==0:
            nngs=4
        else:
            nngs=int(nngs)
    if ndm==None:
        ndm=raw_input("Please enter the number of virtual DMs required: [3] ")
        if len(ndm)==0:
            ndm=3
        else:
            ndm=int(ndm)
    if mpi==None:
        mpi=raw_input("Use MPI?: y/[n] ")
        if len(mpi)==0 or mpi!="y":
            mpi=0
        else:
            mpi=1
    if dirname==None:
        tmp=raw_input("\n\nPlease enter an output directory name:\n[%s] "%dirname)
        if len(tmp)>0:
            dirname=tmp
    if os.path.exists(dirname):
        print("Using existing directory %s"%dirname)
    else:
        os.mkdir(dirname)


    ngs=nngs+nlgs
    txt=""
    for i in range(ngs):
        txt+="%d,"%(i+1)
    if mpi:
        xmltxt=mcaoTxt.replace('groupShare idstr="1" cpu="[]"','groupShare idstr="%s" cpu="[%s]"'%(txt,txt))
    else:#no MPI
        xmltxt=mcaoTxt.replace('groupShare idstr="1"','groupShare idstr="%s"'%txt)


    txt=""
    #remove science
    xmltxt=xmltxt.replace('<groupShare idstr="sci1" cpu="[]" coordlist="[(140, 104), (230, 104), (230, 344), (140, 344)]"/>','')
    #now add the DMs.  Need to:
    #Move tomoRecon, wfscent, science and grouping down.
    #Remove existing DM
    #Break connections
    #Add new DMs
    objlist=xmltxt.split("</simulationObject>")
    newlist=[]
    olddmlist=[]
    gsdmlist=[]
    scidmlist=[]
    for obj in objlist:
        needed=1
        if "304" in obj:#y position of wfscent and science
            obj=obj.replace("304","%d"%(224+80*ndm))
        elif "384" in obj:
            obj=obj.replace("384","%d"%(304+80*ndm))
        if "<groupShare" in obj:
            obj=obj.replace("344","%d"%(264+80*ndm))
        if "tomoRecon" in obj:
            #adjust the connectto...
            txt=""
            for i in range(ndm):
                txt+="%d, "%(i+16)
            obj=obj.replace('connectto="[5, 13]"','connectto="[%s]"'%txt)
            txt=""
            for i in range(ndm):
                txt+="[(129, %d), (129, %d)],\n"%(417+(ndm-1)*80,192+i*80)
            obj=obj.replace("[(129, 417), (129, 192)],\n[(129, 417), (129, 192)],",txt)
        if "iatmos" in obj:
            needed=0
        if "iscrn" in obj:
            obj=obj.split("<simulationObject")[0]#keep the initial bumph.
            obj=obj.replace("doInitialPokeThenRun","doInitialPoke")
        if "wfscent" in obj:
            obj=obj.replace('connectfrom="[5]"','connectfrom="[%d]"'%(16+ndm-1))
        if "science.science" in obj:
            needed=0

        if "science.xinterp_dm" not in obj:
            if needed:
                newlist.append(obj)
        else:
            olddmlist.append(obj)
    dm=olddmlist[0]
    for i in range(ndm):
        txt=dm.replace("80,224","80,%d"%(224+i*80))
        txt=txt.replace('tag="5"','tag="%d"'%(16+i))
        if i==ndm-1:
            connto=7
        else:
            connto=16+i+1
        txt=txt.replace('connectto="[7]"','connectto="[%d]"'%(connto))
        if i==0:
            txt=txt.replace('connectfrom="[3, ','connectfrom="[')
        else:
            connfrom=16+i-1
            txt=txt.replace('connectfrom="[3','connectfrom="[%d'%(connfrom))
        txt=txt.replace("[(129, 417), (129, 192)]","[(129, %d), (129, %d)]"%(417+(ndm-1)*80,192+i*80))
        txt=txt.replace("<sharedTo>\n[13]","<sharedTo>\n[]")
        txt=txt.replace("dm$","dm%dpath$"%i)
        gsdmlist.append(txt)

    gsdmlist[0]=newlist[0]+gsdmlist[0]
    objlist=gsdmlist+newlist[1:]
    xmltxt=string.join(objlist,"</simulationObject>")
    save(xmltxt,dirname,fname+".xml")
    p=util.parseSimXml.parseSimXml(xmltxt)
    pystr=p.makePython()
    save(pystr,dirname,fname+".py")
    if writeParams:
        txt=mcaoParamsTxt.replace("nlgs=4","nlgs=%d"%nlgs)
        txt=txt.replace("nsci=1","nsci=0#Add some science objects if you want to increase the DM field of view beyond that of the guide stars")
        txt=txt.replace("nngs=3","nngs=%d"%nngs)
        txt=txt.replace("ndm=3","ndm=%d"%ndm)
        txt=txt.replace("computeControl=1","computeControl=0\nr.abortAfterPoke=1")
        save(txt,dirname,"paramsPoke.py")
    if mpi:
        execstr="mpirun -np %d python %s.py paramsPoke.py"%(ngs,fname)
        mpitxt="\nThe order of module execution may need editing to maximise performance:\nPlease edit lines starting with 'execOrder' in %s.py to do this.\n(sorry: the logic is difficult to get right for all situations).\n"%fname
    else:
        execstr="python %s.py paramsPoke.py"%fname
        mpitxt=""
    print("""Simulation generated.
To run in %s/, use:
%s
To view the simulation configuration, use:
daspsetup.py %s.xml
To connect to a running simulation (to view WFSs, DMs, etc in action) use:
daspctrl.py
or:
daspanalyse.py

Please edit paramsPoke.py to alter the simulation to meet your requirements.
%s"""%(dirname,execstr,fname,mpitxt))


def makeltao(dirname="ltaoSim",fname="ltao"):
    nlgs=raw_input("Please enter the number of LGS required: [4] ")
    nngs=raw_input("Please enter the number of NGS required: [4] ")
    nsci=raw_input("Please enter the number of science PSFs required: [1] ")
    ndm=raw_input("Please enter the number of virtual DMs required: [3] ")
    if len(nlgs)==0:
        nlgs=4
    else:
        nlgs=int(nlgs)
    if len(nngs)==0:
        nngs=4
    else:
        nngs=int(nngs)
    if len(ndm)==0:
        ndm=3
    else:
        ndm=int(ndm)
    if len(nsci)==0:
        nsci=1
    else:
        nsci=int(nsci)
    mpi=raw_input("Use MPI?: y/[n] ")
    if len(mpi)==0 or mpi!="y":
        mpi=0
    else:
        mpi=1
    tmp=raw_input("\n\nPlease enter an output directory name:\n[%s] "%dirname)
    if len(tmp)>0:
        dirname=tmp
    if os.path.exists(dirname):
        print("Using existing directory %s"%dirname)
    else:
        os.mkdir(dirname)

    #Make the poke simulation:
    makepmx(nlgs,nngs,ndm,mpi,dirname,writeParams=0,fname=fname+"Poke")
    #and a 1 DM MCAO simulation.
    makemcao(nlgs,nngs,nsci,1,mpi,dirname,fname=fname,dopoke=0,addvdm=1)
    save(projectionTxt,dirname,"projection.py")
    ltaoscripttxt="""
#!/bin/sh
python ltaoPoke.py params.py --param="this.tomoRecon.abortAfterPoke=1;this.tomoRecon.reconmxFilename='rmxTomo.fits'" --init="ndm=%d"
python projection.py %d
python ltao.py params.py
"""%(ndm,ndm)
    save(ltaoscripttxt,dirname,"run.sh")
    print("""

****************************************
LTAO Simulation generated.  This does a tomographic reconstruction on %d DMs, and then projects the correction on-axis.

To run in %s/, use (for MPI versions, adjust accordingly):
python ltaoPoke.py params.py --param="this.tomoRecon.abortAfterPoke=1;this.tomoRecon.reconmxFilename='rmxTomo.fits'" --init="ndm=%d"
python projection.py %d
python ltao.py params.py

Or, more simply:
sh run.sh

Some notes about this simulation:
It will first of all take an interaction matrix.
This will then be inverted to give a control matrix, using a simple SVD.
This will then be projected along the on-axis line of sight, to give the LTAO
reconstruction matrix.
This least-squares wavefront reconstructor will then be used for AO correction.

The LGS tip-tilt signal is assumed valid, i.e. is used.  

Several improvements are possible (left as an excercise for the user):
Minimum variance wavefront reconstruction from the interaction matrix.
Ignore the LGS tip-tilt signal, in the user-generated control matrix.
Model-based reconstruction (e.g. FrIM, FEWHA, etc)

Once you have a suitable reconstruction matrix, you do not need to rerun the
poking or projection parts.

    """%(ndm,dirname,ndm,ndm))


def makemoao(dirname="moaoSim",fname="moao"):
    nlgs=raw_input("Please enter the number of LGS required: [4] ")
    nngs=raw_input("Please enter the number of NGS required: [4] ")
    nsci=raw_input("Please enter the number of science PSFs required: [1] ")
    ndm=raw_input("Please enter the number of virtual DMs required: [3] ")
    if len(nlgs)==0:
        nlgs=4
    else:
        nlgs=int(nlgs)
    if len(nngs)==0:
        nngs=4
    else:
        nngs=int(nngs)
    if len(ndm)==0:
        ndm=3
    else:
        ndm=int(ndm)
    if len(nsci)==0:
        nsci=1
    else:
        nsci=int(nsci)
    mpi=raw_input("Use MPI?: y/[n] ")
    if len(mpi)==0 or mpi!="y":
        mpi=0
    else:
        mpi=1
    tmp=raw_input("\n\nPlease enter an output directory name:\n[%s] "%dirname)
    if len(tmp)>0:
        dirname=tmp
    if os.path.exists(dirname):
        print("Using existing directory %s"%dirname)
    else:
        os.mkdir(dirname)
    ngs=nlgs+nngs
    #Make the poke simulation:
    makepmx(nlgs,nngs,ndm,mpi,dirname,writeParams=0,fname=fname+"Poke")
    #make the moao simulation:
    wfstxt=""
    nodetxt=""
    for i in range(ngs):
        wfstxt+="%d,"%(i+1)
        nodetxt+="%d,"%(i+2)
    if mpi==0:
        nodetxt=""
    xmltxt=moaoTxt.replace('<groupShare idstr="1" cpu="[]"','<groupShare idstr="%s" cpu="[%s]"'%(wfstxt,nodetxt))
    scitxt=""
    nodetxt=""
    for i in range(nsci):
        scitxt+="sci%d,"%(i+1)
        nodetxt+="%d,"%(i+2+ngs)
    if mpi==0:
        nodetxt=""
    xmltxt=xmltxt.replace('<groupShare idstr="sci1" cpu="[]"','<groupShare idstr="%s" cpu="[%s]"'%(scitxt,nodetxt))
    save(xmltxt,dirname,fname+".xml")
    
    p=util.parseSimXml.parseSimXml(xmltxt)
    pystr=p.makePython()
    save(pystr,dirname,fname+".py")

    txt=mcaoParamsTxt.replace("nsci=1","nsci=%d"%nsci)
    txt=txt.replace("nlgs=4","nlgs=%d"%nlgs)
    txt=txt.replace("nngs=3","nngs=%d"%nngs)
    txt=txt.replace("ndm=3","ndm=%d"%ndm)
    vdmtxt="""for i in range(nsci):
    #Add the virtual DM projector
    dmInfoList.append(dmInfo('vdmsci%d'%(i+1),["sci%d"%(i+1)],0,nAct,minarea=0.1,actuatorsFrom=["dm%dpath"%x for x in range(ndm)],maxActDist=1.5,decayFactor=0.95,reconLam=lgsLam,closedLoop=0,primaryTheta=atmosGeom.sourceTheta("sci%d"%(i+1)),primaryPhi=atmosGeom.sourcePhi("sci%d"%(i+1))))
    #And the physical MOAO DMs
    dmInfoList.append(dmInfo('dmsci%d'%(i+1),["sci%d"%(i+1)],0,nAct,minarea=0.1,actuatorsFrom="vdmsci%d"%(i+1),maxActDist=1.5,decayFactor=0.,reconLam=lgsLam,closedLoop=0))
"""
    txt=txt.replace("dmOverview=dmOverview(dmInfoList,atmosGeom)",vdmtxt+"dmOverview=dmOverview(dmInfoList,atmosGeom)\nreconIdStr='recon'")
    txt=txt.replace("decayFactor=0.95)","decayFactor=0.,reconLam=lgsLam,closedLoop=0)")
    txt=txt.replace("decayFactor=0.95","decayFactor=0.")
    txt=txt.replace("gainFactor=0.5","gainFactor=1.\nr.abortAfterPoke=1")
    save(txt,dirname,"params.py")
    moaoscripttxt="""
#!/bin/sh
python moaoPoke.py params.py
python moao.py params.py
"""
    save(moaoscripttxt,dirname,"run.sh")
    print("""
To run in %s/, use (for MPI versions, adjust accordingly):
python moaoPoke.py params.py
python moao.py params.py

Or, more simply:
sh run.sh

Some notes about this simulation:
It will first of all take an interaction matrix.
This will then be inverted to give a control matrix, using a simple SVD
(least-squares), which will then be used for wavefront reconstruction.


The LGS tip-tilt signal is assumed valid, i.e. is used.  

Several improvements are possible (left as an excercise for the user):
Minimum variance wavefront reconstruction from the interaction matrix.
Ignore the LGS tip-tilt signal, in the user-generated control matrix.
Model-based reconstruction (e.g. FrIM, FEWHA, etc)

Once you have a suitable reconstruction matrix, you do not need to rerun the
poking simulation.
    """%(dirname))



def makelearn(dirname="learnSim",fname="learn"):
    nlgs=raw_input("Please enter the number of LGS required: [4] ")
    nngs=raw_input("Please enter the number of NGS required: [4] ")
    nsci=raw_input("Please enter the number of science PSFs required: [0] ")
    if len(nlgs)==0:
        nlgs=4
    else:
        nlgs=int(nlgs)
    if len(nngs)==0:
        nngs=4
    else:
        nngs=int(nngs)
    if len(nsci)==0:
        nsci=0
    else:
        nsci=int(nsci)
    mpi=raw_input("Use MPI?: y/[n] ")
    if len(mpi)==0 or mpi!="y":
        mpi=0
    else:
        mpi=1
    tmp=raw_input("\n\nPlease enter an output directory name:\n[%s] "%dirname)
    if len(tmp)>0:
        dirname=tmp
    if os.path.exists(dirname):
        print("Using existing directory %s"%dirname)
    else:
        os.mkdir(dirname)
    ngs=nlgs+nngs

    wfstxt=""
    nodetxt=""
    for i in range(ngs):
        wfstxt+="%d,"%(i+1)
        nodetxt+="%d,"%(i+2)
    if mpi==0:
        nodetxt=""
    xmltxt=learntxt.replace('<groupShare idstr="1" cpu="[]"','<groupShare idstr="%s" cpu="[%s]"'%(wfstxt,nodetxt))
    scitxt=""
    nodetxt=""
    for i in range(nsci):
        scitxt+="sci%d,"%(i+1)
        nodetxt+="%d,"%(i+2+ngs)
    if mpi==0:
        nodetxt=""
    xmltxt=xmltxt.replace('<groupShare idstr="sci1" cpu="[]"','<groupShare idstr="%s" cpu="[%s]"'%(scitxt,nodetxt))
    save(xmltxt,dirname,fname+".xml")
    
    p=util.parseSimXml.parseSimXml(xmltxt)
    pystr=p.makePython()
    save(pystr,dirname,fname+".py")

    txt=mcaoParamsTxt.replace("nsci=1","nsci=%d"%nsci)
    txt=txt.replace("nlgs=4","nlgs=%d"%nlgs)
    txt=txt.replace("nngs=3","nngs=%d"%nngs)
    txt=txt.replace("ndm=3","ndm=0")
    save(txt,dirname,"params.py")
    learnscripttxt="""
#!/bin/sh
python learn.py params.py
"""
    save(learnscripttxt,dirname,"run.sh")
    print("""
To run in %s/, use (for MPI versions, adjust accordingly):
python learn.py params.py

Some notes about this simulation:
It simply records wavefront slope measurements, into files named 
saveOutput*.fits, LGS first, then NGS.  

Uncorrected science PSFs can also be generated, for comparison with other simulations.
    """%(dirname))







projectionTxt="""import sys
import util.dm
util.dm.dmProjectionQuick("params.py",rmx="rmxTomo.fits",rmxOutName="rmx.fits",reconIdStr="recon",initDict={"ndm":int(sys.argv[1])})
"""

##########################################################################

learntxt="""
<aosim>
<simSetup>
<simulationObject cpu="(1, 1)" import="science.iscrn" object="iscrn" pos="80,64" tag="1" shortname="iscrn" pixmap="infScrn.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[3, 11]" connectfrom="[]" textcol="red" idstr="allLayers">
<lines>
[
[],
[],
]
</lines>
<endlines>
[
]
</endlines>
<parentNames>
[]
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.iatmos" object="iatmos" pos="80,144" tag="3" shortname="iatmos" pixmap="infAtmos.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[7]" connectfrom="[1]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['allLayers']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.wfscent" object="wfscent" pos="80,218" tag="7" shortname="wfscent" pixmap="wfscent.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[27]" connectfrom="[3]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.iatmos" object="iatmos" pos="184,144" tag="11" shortname="iatmos" pixmap="infAtmos.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[15]" connectfrom="[1]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['allLayers']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.science" object="science" pos="184,220" tag="15" shortname="science" pixmap="sci.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[]" connectfrom="[11]" textcol="red" idstr="$">
<lines>
[
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<groupShare idstr="1" cpu="[]" coordlist="[(30, 104), (124, 104), (124, 335), (30, 335)]"/>
<groupShare idstr="sci1" cpu="[]" coordlist="[(140, 104), (229, 104), (229, 266), (140, 266)]"/>
<simulationObject cpu="(1, 1)" import="base.saveOutput" object="saveOutput" pos="80,293" tag="27" shortname="save" pixmap="saveOutput.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[]" connectfrom="[7]" textcol="red" idstr="$">
<lines>
[
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
</simSetup>
</aosim>
"""


############################################################################

moaoTxt="""
<aosim>
<simSetup>
<precode>
</precode>
<simulationObject cpu="(1, 1)" import="science.iscrn" object="iscrn" pos="80,64" tag="1" shortname="iscrn" pixmap="infScrn.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[3, 11]" connectfrom="[]" textcol="red" idstr="allLayers">
<lines>
[
[],
[],
]
</lines>
<endlines>
[
]
</endlines>
<parentNames>
[]
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.iatmos" object="iatmos" pos="80,144" tag="3" shortname="iatmos" pixmap="infAtmos.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[7]" connectfrom="[1]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['allLayers']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.wfscent" object="wfscent" pos="80,218" tag="7" shortname="wfscent" pixmap="wfscent.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[9]" connectfrom="[3]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.tomoRecon" object="recon" pos="80,304" tag="9" shortname="tomoRecon" pixmap="xinterp_recon.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[21]" connectfrom="[7]" textcol="red" idstr="recon">
<lines>
[
[(131, 337), (131, 110)],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['$']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.iatmos" object="iatmos" pos="262,142" tag="11" shortname="iatmos" pixmap="infAtmos.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[17]" connectfrom="[1]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['allLayers']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.xinterp_dm" object="dm" pos="262,297" tag="17" shortname="xdm" pixmap="xinterp_dm.xpm" feedback="0" pyname="" groupshare="1" args="" connectto="[15]" connectfrom="[11, 23]" textcol="red" idstr="dm$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
[],
]
</endlines>
<parentNames>
['', '']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.science" object="science" pos="262,372" tag="15" shortname="science" pixmap="sci.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[]" connectfrom="[17]" textcol="red" idstr="$">
<lines>
[
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<groupShare idstr="1" cpu="[]" coordlist="[(30, 104), (121, 104), (121, 266), (30, 266)]"/>
<groupShare idstr="sci1" cpu="[]" coordlist="[(140, 104), (308, 104), (308, 411), (140, 411)]"/>
<simulationObject cpu="(1, 1)" import="science.vdmUser" object="vdmUser" pos="182,142" tag="21" shortname="vdmUser" pixmap="vdmUser.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[23]" connectfrom="[9]" textcol="red" idstr="vdm$">
<lines>
[
[],
]
</lines>
<endlines>
[
[(131, 337), (131, 110)],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="base.fifo" object="fifo" pos="182,220" tag="23" shortname="fifo" pixmap="fifo.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[17]" connectfrom="[21]" textcol="red" idstr="1">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
</simSetup>
</aosim>
"""

#############################################################################


mcaoTxt="""
<aosim>
<simSetup>
<precode>
if not "nopoke" in ctrl.userArgList:ctrl.doInitialPokeThenRun()
</precode>
<simulationObject cpu="(1, 1)" import="science.iscrn" object="iscrn" pos="80,64" tag="1" shortname="iscrn" pixmap="infScrn.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[3, 11]" connectfrom="[]" textcol="red" idstr="allLayers">
<lines>
[
[],
[],
]
</lines>
<endlines>
[
]
</endlines>
<parentNames>
[]
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.iatmos" object="iatmos" pos="80,144" tag="3" shortname="iatmos" pixmap="infAtmos.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[5]" connectfrom="[1]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['allLayers']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.xinterp_dm" object="dm" pos="80,224" tag="5" shortname="xdm" pixmap="xinterp_dm.xpm" feedback="0" pyname="" groupshare="1" args="" connectto="[7]" connectfrom="[3, 9]" textcol="red" idstr="dm$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
[(129, 417), (129, 192)],
]
</endlines>
<parentNames>
['', '']
</parentNames>
<sharedTo>
[13]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.wfscent" object="wfscent" pos="80,304" tag="7" shortname="wfscent" pixmap="wfscent.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[9]" connectfrom="[5]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.tomoRecon" object="recon" pos="80,384" tag="9" shortname="tomoRecon" pixmap="xinterp_recon.xpm" feedback="1" pyname="" groupshare="0" args="" connectto="[5, 13]" connectfrom="[7]" textcol="red" idstr="recon">
<lines>
[
[(129, 417), (129, 192)],
[(129, 417), (129, 192)],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['$']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.iatmos" object="iatmos" pos="180,144" tag="11" shortname="iatmos" pixmap="infAtmos.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[13]" connectfrom="[1]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['allLayers']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.xinterp_dm" object="dm" pos="180,224" tag="13" shortname="xdm" pixmap="xinterp_dm.xpm" feedback="0" pyname="" groupshare="1" args="" connectto="[15]" connectfrom="[9, 11]" textcol="red" idstr="dm$">
<lines>
[
[],
]
</lines>
<endlines>
[
[(129, 417), (129, 192)],
[],
]
</endlines>
<parentNames>
['', '']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.science" object="science" pos="180,304" tag="15" shortname="science" pixmap="sci.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[]" connectfrom="[13]" textcol="red" idstr="$">
<lines>
[
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<groupShare idstr="1" cpu="[]" coordlist="[(30, 104), (120, 104), (120, 344), (30, 344)]"/>
<groupShare idstr="sci1" cpu="[]" coordlist="[(140, 104), (230, 104), (230, 344), (140, 344)]"/>
</simSetup>
</aosim>
"""

#############################################################################


scaoNsciTxt="""<aosim>
<simSetup>
<precode>
if not "nopoke" in ctrl.userArgList:ctrl.doInitialPokeThenRun()
</precode>
<simulationObject cpu="(1, 1)" import="science.iscrn" object="iscrn" pos="64,64" tag="1" shortname="iscrn" pixmap="infScrn.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[3, 4]" connectfrom="[]" textcol="red" idstr="allLayers">
<lines>
[
[],
[],
]
</lines>
<endlines>
[
]
</endlines>
<parentNames>
[]
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.iatmos" object="iatmos" pos="64,144" tag="3" shortname="iatmos" pixmap="infAtmos.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[6]" connectfrom="[1]" textcol="red" idstr="1">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['allLayers']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.iatmos" object="iatmos" pos="164,144" tag="4" shortname="iatmos" pixmap="infAtmos.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[7]" connectfrom="[1]" textcol="red" idstr="$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['allLayers']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.xinterp_dm" object="dm" pos="64,224" tag="6" shortname="xdm" pixmap="xinterp_dm.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[9]" connectfrom="[3, 11]" textcol="red" idstr="dm1">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
[(106, 417), (106, 192)],
]
</endlines>
<parentNames>
['', '']
</parentNames>
<sharedTo>
[7]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.xinterp_dm" object="dm" pos="164,224" tag="7" shortname="xdm" pixmap="xinterp_dm.xpm" feedback="0" pyname="" groupshare="1" args="" connectto="[13]" connectfrom="[4, 11]" textcol="red" idstr="dm$">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
[(122, 417), (122, 192)],
]
</endlines>
<parentNames>
['', '']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.wfscent" object="wfscent" pos="64,303" tag="9" shortname="wfscent" pixmap="wfscent.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[11]" connectfrom="[6]" textcol="red" idstr="1">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.tomoRecon" object="recon" pos="64,384" tag="11" shortname="tomoRecon" pixmap="xinterp_recon.xpm" feedback="1" pyname="" groupshare="0" args="" connectto="[6, 7]" connectfrom="[9]" textcol="red" idstr="recon">
<lines>
[
[(106, 417), (106, 192)],
[(122, 417), (122, 192)],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['1']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.science" object="science" pos="164,304" tag="13" shortname="science" pixmap="sci.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[]" connectfrom="[7]" textcol="red" idstr="$">
<lines>
[
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<groupShare idstr="sci1" cpu="[]" coordlist="[(118, 100), (211, 100), (211, 350), (118, 350)]"/>
</simSetup>
</aosim>
"""

###########################################################################

scao1sciTxt="""
<aosim>
<simSetup>
<precode>
if not "nopoke" in ctrl.userArgList:ctrl.doInitialPokeThenRun()
</precode>
<simulationObject cpu="(1, 1)" import="science.iscrn" object="iscrn" pos="70,64" tag="1" shortname="iscrn" pixmap="infScrn.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[4]" connectfrom="[]" textcol="red" idstr="allLayers">
<lines>
[
[],
]
</lines>
<endlines>
[
]
</endlines>
<parentNames>
[]
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.iatmos" object="iatmos" pos="70,144" tag="4" shortname="iatmos" pixmap="infAtmos.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[6]" connectfrom="[1]" textcol="red" idstr="1">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['allLayers']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.xinterp_dm" object="dm" pos="70,224" tag="6" shortname="xdm" pixmap="xinterp_dm.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[8, 12]" connectfrom="[4, 10]" textcol="red" idstr="dm1">
<lines>
[
[],
[],
]
</lines>
<endlines>
[
[],
[(30, 424), (30, 184)],
]
</endlines>
<parentNames>
['', '']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.wfscent" object="wfscent" pos="70,304" tag="8" shortname="wfscent" pixmap="wfscent.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[10]" connectfrom="[6]" textcol="red" idstr="1">
<lines>
[
[],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.tomoRecon" object="recon" pos="70,384" tag="10" shortname="tomoRecon" pixmap="xinterp_recon.xpm" feedback="1" pyname="" groupshare="0" args="" connectto="[6]" connectfrom="[8]" textcol="red" idstr="recon">
<lines>
[
[(30, 424), (30, 184)],
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['1']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
<simulationObject cpu="(1, 1)" import="science.science" object="science" pos="170,304" tag="12" shortname="science" pixmap="sci.xpm" feedback="0" pyname="" groupshare="0" args="" connectto="[]" connectfrom="[6]" textcol="red" idstr="sci1">
<lines>
[
]
</lines>
<endlines>
[
[],
]
</endlines>
<parentNames>
['']
</parentNames>
<sharedTo>
[]
</sharedTo>
</simulationObject>
</simSetup>
</aosim>
"""

###########################################################################

mcaoParamsTxt="""
wfs_nsubx=6 #Number of subaps
tstep=1/250.#Simulation timestep in seconds (250Hz).
AOExpTime=40.#40 seconds exposure (use --iterations=xxx to modify)
phasesize=8#number of phase pixels per sub-aperture.
npup=wfs_nsubx*phasesize#Number of phase points across the pupil
telDiam=4.2*wfs_nsubx/6.#Telescope diameter
telSec=telDiam/7.#Central obscuration
ntel=npup#Telescope diameter in pixels
nAct=wfs_nsubx+1#Number of actuators across the DM
ngsLam=640.#NGS wavelength
lgsLam=589.#LGS wavelength
sciLam=1650.#Science wavelength
lgsAsterismRadius=60.#arcseconds
ngsAsterismRadius=90.#arcseconds
nsci=1
nlgs=4
nngs=3
if hasattr(this.globals,"ndm"):
    ndm=this.globals.ndm
else:
    ndm=3
import util.tel
#Create a pupil function
pupil=util.tel.Pupil(npup,ntel/2,ntel/2*telSec/telDiam)

#Create the WFS overview
lgsalt=90000.#90km sodium layer
import util.guideStar
import util.elong
#Create the LGS PSFs (elongated).  There are many ways to do this - this is a simple one.
psf=util.elong.make(spotsize=phasesize*4,nsubx=wfs_nsubx,wfs_n=phasesize,lam=lgsLam,telDiam=telDiam,telSec=telSec,beacon_alt=lgsalt,beacon_depth=10000.,launchDist=0.,launchTheta=0.,pup=pupil)[0]

sourceList=[]
wfsDict={}
for i in range(nlgs):#60 arcsec off-axis
    id="%d"%(i+1)
    wfsDict[id]=util.guideStar.LGS(id,wfs_nsubx,lgsAsterismRadius,i*360./nlgs,lgsalt,phasesize=phasesize,minarea=0.5,sig=1e6,sourcelam=lgsLam,reconList=["recon"],pupil=pupil,launchDist=0,launchTheta=0,lgsPsf=psf)
    sourceList.append(wfsDict[id])
for i in range(nngs):#90 arcsec off-axis
    id="%d"%(i+1+nlgs)
    wfsDict[id]=util.guideStar.NGS(id,wfs_nsubx,ngsAsterismRadius,i*360./nngs,phasesize=phasesize,minarea=0.5,sig=1e6,sourcelam=ngsLam,reconList=["recon"],pupil=pupil)
    sourceList.append(wfsDict[id])
wfsOverview=util.guideStar.wfsOverview(wfsDict)

#Create a Science overview.
import util.sci
sciDict={}
for i in range(nsci):
    id="sci%d"%(i+1)
    sciDict[id]=util.sci.sciInfo(id,i*10.,0.,pupil,sciLam,phslam=sciLam)
    sourceList.append(sciDict[id])
sciOverview=util.sci.sciOverview(sciDict)
#Create the atmosphere object and source directions.
from util.atmos import geom,layer,source
atmosDict={}
nlayer=10 #10 atmospheric layer
layerList={"allLayers":["L%d"%x for x in range(nlayer)]}
strList=[0.5]+[0.5/(nlayer-1.)]*(nlayer-1)#relative strength of the layers
hList=range(0,nlayer*1000,1000)#height of the layers
vList=[10.]*nlayer#velocity of the layers
dirList=range(0,nlayer*10,10)#direction (degrees) of the layers
for i in range(nlayer):
 atmosDict["L%d"%i]=layer(hList[i],dirList[i],vList[i],strList[i],10+i)

l0=10. #outer scale
r0=0.137 #fried's parameter
atmosGeom=geom(atmosDict,sourceList,ntel,npup,telDiam,r0,l0)


#Create the DM object.
from util.dm import dmOverview,dmInfo
import numpy
if ndm>1:
    dmHeight=numpy.arange(ndm)*(hList[-1]/(ndm-1.))
else:
    dmHeight=[0]
dmInfoList=[]
for i in range(ndm):
    dmInfoList.append(dmInfo('dm%dpath'%i,[x.idstr for x in sourceList],dmHeight[i],nAct,minarea=0.1,actuatorsFrom="recon",pokeSpacing=(None if wfs_nsubx<20 else 10),maxActDist=1.5,decayFactor=0.95))
dmOverview=dmOverview(dmInfoList,atmosGeom)

#reconstructor
this.tomoRecon=new()
r=this.tomoRecon
r.rcond=0.05#condtioning value for SVD
r.recontype="pinv"#reconstruction type
r.pokeval=1.#strength of poke
r.gainFactor=0.5#Loop gain
r.computeControl=1#To compute the control matrix after poking
r.reconmxFilename="rmx.fits"#control matrix name (will be created)
r.pmxFilename="pmx.fits"#interation matrix name (will be created)
"""


##############################################################################

scaoParamsTxt="""
wfs_nsubx=6 #Number of subaps
tstep=1/250.#Simulation timestep in seconds (250Hz).
AOExpTime=40.#40 seconds exposure (use --iterations=xxx to modify)
npup=wfs_nsubx*8#Number of phase points across the pupil
telDiam=4.2*wfs_nsubx/6.#Telescope diameter
telSec=telDiam/7.#Central obscuration
ntel=npup#Telescope diameter in pixels
nAct=wfs_nsubx+1#Number of actuators across the DM
ngsLam=640.#NGS wavelength
sciLam=1650.#Science wavelength
nsci=1
import util.tel
#Create a pupil function
pupil=util.tel.Pupil(npup,ntel/2,ntel/2*telSec/telDiam)

#Create the WFS overview
import util.guideStar
wfsDict={"1":util.guideStar.NGS("1",wfs_nsubx,0.,0.,phasesize=npup/wfs_nsubx,\
                                minarea=0.5,sig=1e6,sourcelam=ngsLam,\
                                reconList=["recon"],pupil=pupil)}
wfsOverview=util.guideStar.wfsOverview(wfsDict)

#Create a Science overview.
import util.sci
sciDict={}
if nsci==1:
 phslam=ngsLam
else:
 phslam=sciLam
for i in range(nsci):
 sciDict["sci%d"%(i+1)]=util.sci.sciInfo("sci%d"%(i+1),i*10.,0.,pupil,sciLam,phslam=phslam)
 sciOverview=util.sci.sciOverview(sciDict)

#Create the atmosphere object and source directions.
from util.atmos import geom,layer,source
atmosDict={}
nlayer=2 #2 atmospheric layers
layerList={"allLayers":["L%d"%x for x in range(nlayer)]}
strList=[0.9]+[0.1]*(nlayer-1)#relative strength of the layers
hList=range(0,nlayer*1000,1000)#height of the layers
vList=[10.]*nlayer#velocity of the layers
dirList=[0.]*nlayer#direction (degrees) of the layers
for i in range(nlayer):
 atmosDict["L%d"%i]=layer(hList[i],dirList[i],vList[i],strList[i],10+i)
sourceList=[]
#the wfs
sourceList.append(wfsOverview.getWfsByID("1"))

#and psf
for i in range(nsci):
 sourceList.append(sciOverview.getSciByID("sci%d"%(i+1)))
l0=10. #outer scale
r0=0.137 #fried's parameter
atmosGeom=geom(atmosDict,sourceList,ntel,npup,telDiam,r0,l0)


#Create the DM object.
from util.dm import dmOverview,dmInfo
dmInfoList=[dmInfo('dm',[x.idstr for x in sourceList],0.,nAct,minarea=0.1,actuatorsFrom="recon",\
                   pokeSpacing=(None if wfs_nsubx<20 else 10),maxActDist=1.5,decayFactor=0.95)]
dmOverview=dmOverview(dmInfoList,atmosGeom)

#reconstructor
this.tomoRecon=new()
r=this.tomoRecon
r.rcond=0.05#condtioning value for SVD
r.recontype="pinv"#reconstruction type
r.pokeval=1.#strength of poke
r.gainFactor=0.5#Loop gain
r.computeControl=1#To compute the control matrix after poking
r.reconmxFilename="rmx.fits"#control matrix name (will be created)
r.pmxFilename="pmx.fits"#interation matrix name (will be created)
"""



if __name__=="__main__":
    makesim()
