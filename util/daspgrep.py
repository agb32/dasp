#!/usr/bin/env python
import sys
import os
import string

paramList=["strehl","strehlpeak","inbox","d50","fwhm","rms","rmserr","rmserrerr"]

def getArgs(args):
    glist=[]#list of strings required for selection
    flist=[]#list of files
    ilist=[]#outputs - strehl, etc.
    exclude=[]#list of strings to exclude
    vallist=[]#values to be printed
    valcharlist=[]#tuples of (N,text) where N is the number of characters after finding text to print.
    printid=0
    printdict=0
    printall=0
    printfile=0
    printindex=0
    precision=0
    space=0
    printval=0
    for a in args:
        if a[:7]=="--grep=":
            glist.append(a[7:])
        elif a[:6]=="--help":
            print "Usage: --grep=STRING --file=FILENAME --param=PARAMETER [--printid --printidBefore --printall --printdict --printall --printfile --printindex --space --precision --exclude=xxx --val=xxx --printval"
            print "Or:  grep string filename parameter"
            print "Note, strehl and inbox can be prefixed with % to return in percentage, eg %strehl %inbox0.1"
            print "Exclude parameter is a text string to be excluded"
            print "val parameter is a text string after which the next value is printed"
            sys.exit(0)
        elif a[:7]=="--file=":
            flist.append(a[7:])
        elif a[:8]=="--param=":
            ilist.append(a[8:])
        elif a[:9]=="--printid":
            if a[:15]=="--printidBefore":
                printid=-1
            else:
                printid=1
            if a[:11]=="--printid=0":
                printid=0
        elif a[:11]=="--printdict":
            printdict=1
            if a[:13]=="--printdict=0":
                printdict=0
        elif a[:10]=="--printall":
            printall=1
            if a[:12]=="--printall=0":
                printall=0
        elif a[:11]=="--printfile":
            printfile=1
            if a[:13]=="--printfile=0":
                printfile=0
        elif a[:12]=="--printindex":
            printindex=1
            if a[:14]=="--printindex=0":
                printindex=0
        elif a[:7]=="--space":
            space=1
            if a=="--space=0":
                space=0
        elif a[:11]=="--precision":
            precision=1
        elif a[:10]=="--exclude=":
            exclude.append(a[10:])
        elif a[:6]=="--val=":
            vallist.append(a[6:])
        elif a[:5]=="--val":
            valcharlist.append((int(a[5:a.index("=")]),a[a.index("=")+1:]))
        elif a[:10]=="--printval":
            printval=1
        else:
            if os.path.exists(a):#is it a filename?
                flist.append(a)
            else:
                got=0
                for p in paramList:
                    if p==string.lower(a)[:len(p)] or (a[0]=="%" and p==string.lower(a[1:])[:len(p)]):#its a parameter
                        ilist.append(a)
                        got=1
                if got==0:#its a grep string
                    glist.append(a)
    return glist,flist,ilist,printid,printdict,printall,printfile,printindex,space,precision,exclude,vallist,printval,valcharlist


def grep(glist,flist,ilist,printid=0,printdict=0,printall=0,printfile=0,printindex=0,space=0,precision=0,exclude=[],vallist=[],printval=0,valcharlist=[]):
    outtxt=""
    pretxt=""
    cnt=0
    fillchr="\t"
    if space:
        fillchr=" "
    for f in flist:
        lines=open(f).readlines()
        for line in lines:
            ok=1
            for g in glist:
                if g not in line:
                    ok=0
                    break
            for e in exclude:
                if e in line:
                    ok=0
                    break
            if ok:#select the parameters now...
                try:
                    indx=line.index("{")
                    dicttxt=line[indx:]
                except:
                    dicttxt="{}"
                try:
                    indx1=dicttxt.index("}")
                    dicttxt=dicttxt[:indx1+1]
                except:
                    dicttxt="{}"
                try:
                    sciDict=eval(dicttxt)
                except:
                    sciDict={}
                try:
                    indx1=line.index("RMS: ")
                    rmsList=line[indx1+5:].split()
                    sciDict["RMS"]=float(rmsList[0])
                    sciDict["RMSErr"]=float(rmsList[2])
                    sciDict["RMSErrErr"]=float(rmsList[4])
                except:
                    pass
                if printall:
                    ilist=sciDict.keys()
                txt=""
                if printindex:
                    txt+="%s%d"%(fillchr,cnt)
                cnt+=1
                if printfile:
                    txt+="%s%s"%(fillchr,f)
                if printid==1:
                    txt+="%s%s"%(fillchr,line[:indx])
                elif printid==-1:
                    pretxt+="%s%s\n"%(fillchr,line[:indx])
                for v in vallist:
                    try:
                        indx=line.index(v)+len(v)
                    except:
                        indx=None
                    if not printval:
                        v=""
                    if indx!=None:
                        pos=1
                        val=None
                        while 1:
                            try:
                                val=eval(line[indx:indx+pos])
                                pos+=1
                            except:
                                break
                        txt+="%s%s%s"%(fillchr,v,str(val))
                    else:
                        txt+="%s%sNone"%(fillchr,v)
                for n,v in valcharlist:
                    try:
                        indx=line.index(v)+len(v)
                    except:
                        indx=None
                    if not printval:
                        v=""
                    if indx!=None:
                        txt+="%s%s%s"%(fillchr,v,line[indx:indx+n])
                    else:
                        txt+="%s%sNONE"%(fillchr,v)
                for param in ilist:
                    if param[0]=="%":
                        m=100.
                        key=param[1:]
                    else:
                        key=param
                        m=1.
                    #if param in ["strehl","strehlPeak"] or param[:5]=="inbox":
                    #    m=100.
                    #else:
                    #    m=1.
                    if sciDict.has_key(key):
                        if printdict:
                            txt+="%s%s"%(fillchr,key)
                        try:
                            if precision==0:
                                txt+="%s%.3g"%(fillchr,sciDict[key]*m)
                            else:
                                txt+="%s%g"%(fillchr,sciDict[key]*m)
                        except:
                            txt+="%s%s"%(fillchr,str(sciDict[key]))

                outtxt+="%s\n"%txt[len(fillchr):]
    return pretxt+outtxt


if __name__=="__main__":
    glist,flist,ilist,printid,printdict,printall,printfile,printindex,space,precision,exclude,vallist,printval,valcharlist=getArgs(sys.argv[1:])
    txt=grep(glist,flist,ilist,printid,printdict,printall,printfile,printindex,space,precision,exclude,vallist,printval,valcharlist)
    print txt
