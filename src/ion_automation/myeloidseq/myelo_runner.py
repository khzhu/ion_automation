import sys
from ion_automation.myeloidseq.myelo_worker import myeloseq

if len(sys.argv) != 2:
    print("Worksheet was not provided, try again!")
else:
    myelo_runner = myeloseq("/home/ionadmin/ion_config.conf")
    myelo_runner.workbook = sys.argv[1]
    myelo_runner.start()