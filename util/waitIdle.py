import subprocess
import time
import sys
#uptime=subprocess.Popen("uptime",stdout=subprocess.PIPE).stdout.read()
def run(ignoremem=0):
    if ignoremem==0:
        print "Waiting for load to be less than 0.1, and more than half of memory free"
    else:
        print "Waiting for load to be less than 0.1"
    okay=0
    while okay==0:
        okay=1
        txt=time.strftime("%H:%M")+": load: "
        top=subprocess.Popen(["top","-n 1"],stdout=subprocess.PIPE).stdout.readlines()
        load=top[0][top[0].index("load average:"):].split()[2:5]
        loadlist=[]
        ps=float(subprocess.Popen("ps -eo pmem,pcpu,ucmd --no-header | sort -r -n | head -n 1",stdout=subprocess.PIPE,shell=True).stdout.read().strip().split()[0])
        for l in load:
            i=0
            while l[i] in ["0","1","2","3","4","5","6","7","8","9","."]:
                i+=1
            loadlist.append(float(l[:i]))
            if loadlist[-1]>0.1:
                okay=0
            txt+="%.2f, "%loadlist[-1]
        ind=2
        tt=top[3].split()
        if tt[0]=="KiB":
            ind+=1
        mem=top[3].split()[ind::2]#total, used, free, buffers
        fracmem=float(mem[1][:-1])/float(mem[0][:-1])
        #if fracmem>0.5:
        #    okay=0
        if ps>1 and ignoremem==0:#more than 1% of memory used by a process - so its probably a simulation...
            okay=0
        if ignoremem==0:
            txt+="used mem: %.1f%% (max per process %.1f%%)   \r"%(fracmem*100,ps)
        else:
            txt+="      \r"
        #print "load: %s used mem: %s\r"%(loadlist,fracmem),
        print txt,
        sys.stdout.flush()
        if okay==0:
            time.sleep(60)

if __name__=="__main__":
    ignoremem=0
    if "--ignore-mem" in sys.argv:
        ignoremem=1
    run(ignoremem)
