#! /bin/bash 
import os, sys, time
from testNameList import nameList, coreNumList
# This script will perform tests on default test cases, listed in testNameList.

os.system('rm -rf test')
os.system('rm -rf bin/eqquasi')
os.system('mkdir test')

os.system('./install.eqquasi.sh -m ubuntu')
os.chdir('test')

startTime = time.time()
def runTest(testDir, compSet, coreNum):
    cmd = 'create.newcase '+testDir+' '+compSet
    os.system(cmd)
    os.chdir(testDir)
    os.system('./case.setup')
    os.system('bash run.sh')
    #os.system('python3 plotRuptureDynamics')
    os.chdir('..')
    
for testName, coreNum in zip(nameList, coreNumList):
    runTest(testName, testName, coreNum)

os.chdir('..')
rc = os.system('python3 check.test.py')

endTime = time.time()

print('Time consumed for all the tests are ', endTime-startTime, ' s')

# Propagate check.test.py's verdict so CI can fail on a regression.
sys.exit(1 if os.waitstatus_to_exitcode(rc) else 0)
