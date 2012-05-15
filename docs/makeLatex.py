"""Python code to make the simulation documentation in latex format.  This can then be converted to PDF or postscript using standard unix commands (see the Makefile)."""
import os,stat,sys,os.path,string


def scandir(prefix):
    pylist=[]
    dirs=os.listdir(prefix)
    def _testdir(name,path):
       for x in path:
          try:
             if os.path.samefile(name,x):
                return 1
          except OSError:
            pass
       return 0
    if "__init__.py" in dirs or _testdir(prefix, sys.path):
        #If its in the python path, or has an __init__.py file present (so that it can be obtained from the python path) do this...
        for d in dirs:
            s=os.stat(prefix+d)
            if stat.S_ISDIR(s.st_mode):
                pylist+=scandir(prefix+d+"/")
            else:
                if d[-3:]==".py" and d!="__init__.py":
                    module=prefix+d
                    module=module.replace("/",".")[:-3]
                    module=module.replace("...","")
                    
                    pylist.append(module)
    return pylist



if __name__=="__main__":
    nolatex=os.system("latex -v >& /dev/null")
    if nolatex!=0:#no latex present
        print "Latex needed to generate documentation.  Present on node 6 of the Cray"
    else:
        modlist=scandir("../")
        removeList=["util.pyfits","docs.makeLatex","base.joiner2","base.splitter"]
        for rem in removeList:
            if rem in modlist:
                modlist.remove(rem)
        modstr=string.join(modlist)
        print "Documenting modules",modstr
        if os.system("epydoc --latex -o modules %s"%modstr)!=0:
            print "Error: epydoc not found or no modules specified"
            raise Exception("Epydoc not found or no modules specified")
        simlines=open("simapi-orig.tex").readlines()
                
        inclines=open("modules/api.tex").readlines()
        f=open("simapi.tex","w")

        includes=[]
        for line in inclines:
            if line[:9]=="\\include{":
                includes.append(line[:9]+"modules/"+line[9:])
##             if line=="\\begin{document}\n":
##                 f.write(line)

##                 f.write("\\setlength{\\headsep}{2cm}\n\\setlength{\\voffset}{-1cm}\n")
##             elif line=="\\renewcommand{\\sectionmark}[1]{\\markboth{#1}{}}\n":
##                 f.write("\\renewcommand{\\sectionmark}[1]{\\markboth{#1}{AO simulation documentation}}\n")
##             elif line=="\\renewcommand{\\subsectionmark}[1]{\\markright{#1}}\n":
##                 f.write("\\renewcommand{\\subsectionmark}[1]{\\markright{AO simulation documentation}}\n")
##             elif line=="\\title{}\n":
##                 f.write("\\title{AO simulation documentation}\n")
##             elif line=="\\author{API Documentation}\n":
##                 f.write("\\author{A. G. Basden}\n")
##             else:
##                 f.write(line)
        

        for line in simlines:
            f.write(line)
            if line=="%Add includes here\n":
                for inc in includes:
                    f.write(inc)
        f.close()


        print "#"*79
        print "Making API documentation"
        print "#"*79
        if len(sys.argv)>1:
            os.system("make simapi.pdf")

