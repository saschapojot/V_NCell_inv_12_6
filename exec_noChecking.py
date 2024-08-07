import subprocess
from decimal import Decimal


#this scrip executes mc without checking statistics
def format_using_decimal(value):
    # Convert the float to a Decimal
    decimal_value = Decimal(value)
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

T=1
unitCellNum=3


TStr=format_using_decimal(T)
#############################################
#launch mc, i.e., giving initial conditions

launchResult=subprocess.run(["python3", "launch_one_run.py", "./dataAllUnitCell"+str(unitCellNum)+"/row0/T"+TStr+"/run_T"+str(T)+".mc.conf"])
print(launchResult.stdout)
if launchResult.returncode!=0:
    print("error in launch one run: "+str(launchResult.returncode))
#############################################


#############################################
#make mc
targetName="run_mc"
compileErrCode=10
cmake_process=subprocess.Popen(["cmake","."], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
while True:
    output = cmake_process.stdout.readline()
    if output == '' and cmake_process.poll() is not None:
        break
    if output:
        print(output.strip())
stdout, stderr = cmake_process.communicate()
if stdout:
    print(stdout.strip())
if stderr:
    print(stderr.strip())



#############################################

#############################################
#run executable
cppExecutable="./run_mc"
process = subprocess.Popen([cppExecutable, "./dataAllUnitCell"+str(unitCellNum)+"/row0/T"+TStr+"/cppIn.txt" ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
while True:
    output = process.stdout.readline()
    if output == '' and process.poll() is not None:
        break
    if output:
        print(output.strip())
stdout, stderr = process.communicate()
if stdout:
    print(stdout.strip())
if stderr:
    print(stderr.strip())

#############################################
