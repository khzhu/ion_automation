import sys
from ion_automation.oncomine_solid.ion_worker import oncomine_solid

if len(sys.argv) != 2:
    print("Worksheet was not provided, try again!")
else:
    ion_worker = oncomine_solid("/home/ionadmin/ion_config.conf")
    ion_worker.workbook = sys.argv[1]
    ion_worker.run()