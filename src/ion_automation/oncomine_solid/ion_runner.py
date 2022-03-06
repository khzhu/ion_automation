from ion_automation.oncomine_solid.ion_worker import oncomine_solid

ion_worker = oncomine_solid("/home/ionadmin/ion_config.conf")
ion_worker.workbook = "/mnt/Z_drive/Molecular/IonTorrent/oncosolid_autoreport/worksheet/22-MGON3.xlsm"
ion_worker.run()