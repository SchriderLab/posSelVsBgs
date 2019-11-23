import sys, os

"""
Simple helper script for using SLURM to launch an array job
"""

def runCmdAsJobArrayWithoutWaitingWithLog(cmd, jobName, launchFile, wallTime, qName, mbMem, logFile, repRange):
    with open(launchFile,"w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=%s\n" %(jobName))
        f.write("#SBATCH --time=%s\n" %(wallTime))
        f.write("#SBATCH --partition=%s\n" %(qName))
        f.write("#SBATCH --output=%s\n" %(logFile))
        f.write("#SBATCH --mem=%s\n" %(mbMem))
        f.write("#SBATCH --requeue\n")
        f.write("#SBATCH --export=ALL\n")
        f.write("\n%s\n" %(cmd))
    os.system("sbatch --array=%s %s" %(repRange, launchFile))

def main():
    cmd, jobName, launchFile, wallTime, qName, mbMem, logFile, repRange = sys.argv[1:]
    runCmdAsJobArrayWithoutWaitingWithLog(cmd, jobName, launchFile, wallTime, qName, mbMem, logFile, repRange)

if __name__ == "__main__":
    main()
